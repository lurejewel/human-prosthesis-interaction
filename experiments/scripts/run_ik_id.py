"""Automated OpenSim IK + ID pipeline with 6 Hz GRF pre-filtering.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\run_ik_id.py --masked"


Steps:
  1. Filter ALL 18 GRF data columns (force / COP / torque for both belts) with a
     4th-order zero-phase Butterworth low-pass at 6 Hz  ->  level_walking_grf_6Hz.mot
     + level_walking_grf_6Hz.xml (external loads pointing to the filtered file).
  2. Inverse Kinematics (Setup_IK.xml)     -> overwrite level_walking_ik.mot
  3. Inverse Dynamics   (Setup_InverseDynamics.xml)
     - uses the FILTERED GRF (level_walking_grf_6Hz.xml)
     - keeps the 6 Hz coordinate filtering and forces_to_exclude=Muscles from the setup
     -> overwrite ResultsInverseDynamics/level_walking_id.sto

Optional flags: --skip-filter, --skip-ik, --skip-id, --masked
    --masked : replace the plain 6 Hz filter with a physics-constrained one
               (per-belt stance mask -> Tukey ramps -> vy>=0, cross-terms zeroed
               during swing, raw COP kept in swing). The ID step then uses the
               masked external loads (level_walking_grf_masked.xml).

File map (everything under experiments/data/SQR_walking/):
  Inputs (must exist, checked by check_inputs):
    level_walking_grf.mot    raw GRF; time + 18 cols: belt1 (1-9), belt2 (10-18),
                             each belt = vx vy vz | px py pz | tx ty tz
    level_walking_grf.xml    ExternalLoads template; only its <datafile> is
                             repointed (regex) to the filtered .mot
    SQR_simbody.osim         model used by both IK and ID
    level_walking.trc        marker trajectories for IK
    Setup_IK.xml / Setup_InverseDynamics.xml   tool setups; ID keeps its 6 Hz
                             coordinate filter and forces_to_exclude=Muscles
  Outputs:
    level_walking_grf_6Hz.mot / .xml          filtered GRF + repointed loads
    level_walking_grf_masked.mot / .xml       masked variant (only with --masked)
    level_walking_grf_*_comparison_*.png      raw-vs-filtered / raw-vs-masked plots
    level_walking_ik.mot                      IK result (previous copy -> .bak)
    ResultsInverseDynamics/level_walking_id.sto  ID result (previous copy -> .bak)

Key parameters (top of file):
    GRF_FILTER_FC = 6.0 Hz   matches the 6 Hz coordinate filter in Setup_InverseDynamics.xml
    TIME_RANGE    = 180-240 s  analysis window; raw files are longer recordings
    MASKED_*      contact/debounce/taper thresholds for the --masked mode

Behaviour notes:
    - The script chdir()s into DATA_DIR so OpenSim resolves relative paths the
      same way the GUI does; that is why outputs overwrite files in place.
    - Sampling rate is estimated from the median time step (fs = 1/median(dt)).
    - run_id() always uses the FILTERED external loads (6 Hz, or masked with
      --masked), never the raw GRF.
    - Dependencies: numpy, scipy, opensim (required); matplotlib only for the
      comparison plots (skipped silently if missing).
"""
import argparse
import glob
import os
import re
import shutil
import sys

import numpy as np
from scipy.ndimage import label
from scipy.signal import butter, filtfilt

import opensim as osim

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

GRF_RAW_MOT = os.path.join(DATA_DIR, "level_walking_grf.mot")
GRF_RAW_XML = os.path.join(DATA_DIR, "level_walking_grf.xml")
GRF_FILT_MOT = os.path.join(DATA_DIR, "level_walking_grf_6Hz.mot")
GRF_FILT_XML = os.path.join(DATA_DIR, "level_walking_grf_6Hz.xml")

MODEL  = os.path.join(DATA_DIR, "SQR_simbody.osim")
TRC    = os.path.join(DATA_DIR, "level_walking.trc")

IK_SETUP = os.path.join(DATA_DIR, "Setup_IK.xml")
IK_OUT   = os.path.join(DATA_DIR, "level_walking_ik.mot")

ID_SETUP    = os.path.join(DATA_DIR, "Setup_InverseDynamics.xml")
ID_RESULTS  = os.path.join(DATA_DIR, "ResultsInverseDynamics")
ID_OUT      = os.path.join(ID_RESULTS, "level_walking_id.sto")
GRF_FILTER_FC = 6.0     # Hz, matches the 6 Hz coordinate filter in Setup_InverseDynamics.xml
TIME_RANGE    = (180.0, 240.0)

# Masked filtering (physical-constraint) configuration
MASKED_THR          = 10.0   # N, vertical-force contact threshold
MASKED_MIN_STANCE_S = 0.10   # drop stance segments shorter than this
MASKED_MIN_GAP_S    = 0.04   # fill swing gaps shorter than this
MASKED_TAPER_MS     = 25.0   # per-end cosine ramp (50 ms total); >=1 kinematics frame (100 Hz)
GRF_MASKED_MOT = os.path.join(DATA_DIR, "level_walking_grf_masked.mot")
GRF_MASKED_XML = os.path.join(DATA_DIR, "level_walking_grf_masked.xml")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def read_table(path):
    """Read an OpenSim .mot/.sto file.

    Returns (header_lines, column_labels, data_ndarray). The column-labels
    line sits right after 'endheader' in the OpenSim text-table format.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    header = []
    i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("endheader"):
        header.append(lines[i].rstrip("\r\n"))
        i += 1
    if i < len(lines):
        header.append("endheader")
        i += 1

    # Skip the column-labels line if present (its first token is not numeric)
    col_labels = None
    while i < len(lines):
        ln = lines[i]
        if ln.strip():
            try:
                float(ln.lstrip().split()[0])
            except ValueError:
                col_labels = ln.strip()
                i += 1
            break
        i += 1

    data_lines = [ln for ln in lines[i:] if ln.strip()]
    data = np.loadtxt(data_lines)
    return header, col_labels, data


def write_table(path, header, col_labels, data):
    """Write an OpenSim .mot/.sto file preserving the original header format."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        for h in header:
            f.write(h + "\n")
        if col_labels:
            f.write(col_labels + "\n")
        for row in data:
            f.write("\t".join("{:.8f}".format(v) for v in row) + "\n")


def column_names(col_labels):
    """Split the column-labels line into a list of names."""
    if not col_labels:
        return []
    return [c.strip() for c in col_labels.split()]


def backup_if_exists(path):
    if os.path.exists(path):
        bak = path + ".bak"
        shutil.copy2(path, bak)
        print("   backed up existing file -> {}".format(bak))


def check_inputs():
    required = [GRF_RAW_MOT, GRF_RAW_XML, MODEL, TRC, IK_SETUP, ID_SETUP]
    missing = [p for p in required if not os.path.exists(p)]
    if missing:
        sys.exit("Missing input files:\n  " + "\n  ".join(missing))
    print("All input files present.")


# 9 GRF components per belt, in file column order:
#   col 1-9  = belt1 vx vy vz | px py pz | tx ty tz
#   col 10-18= belt2 vx vy vz | px py pz | tx ty tz
GRF_COMPONENTS = ["vx (N)", "vy (N)", "vz (N)",
                  "px (m)", "py (m)", "pz (m)",
                  "tx (Nm)", "ty (Nm)", "tz (Nm)"]


def plot_grf_comparison(t, raw, filt, out_path, tlim=None, downsample=1,
                        filt_label="belt2 filt", filt_desc=None):
    """Overlay raw vs filtered GRF, 3x3 subplots (one per component).

    Only belt 2 (ground_force2_v / ground_force2_p / ground_torque2_,
    file columns 10-18) is plotted: each subplot shows 2 curves (raw, filtered).
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("   [skip plot] matplotlib not available: " + out_path)
        return

    sel = np.ones(len(t), dtype=bool)
    if tlim is not None:
        sel = (t >= tlim[0]) & (t <= tlim[1])
    tt = t[sel][::downsample]
    raw_s = raw[sel][::downsample]
    filt_s = filt[sel][::downsample]

    fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharex=True)
    for k, comp in enumerate(GRF_COMPONENTS):
        ax = axes[k // 3][k % 3]
        c2 = 10 + k                      # belt2 (ground_force2_*) column offset
        ax.plot(tt, raw_s[:, c2],  color="C0", lw=0.7, label="belt2 raw")
        ax.plot(tt, filt_s[:, c2], color="C1", lw=0.7, label=filt_label)
        ax.set_title(comp)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)
    if filt_desc is None:
        filt_desc = "{} Hz filtered".format(int(GRF_FILTER_FC))
    fig.suptitle("GRF raw vs {} ({})".format(filt_desc, os.path.basename(out_path)))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print("   saved plot -> {}".format(out_path))


# ---------------------------------------------------------------------------
# Masked-filtering helpers
# ---------------------------------------------------------------------------
def debounce_mask(on, fs, min_stance_s=0.10, min_gap_s=0.04):
    """Clean a boolean contact mask: drop stance segments shorter than
    `min_stance_s` and fill swing gaps shorter than `min_gap_s`."""
    on = on.copy()
    min_stance_n = max(int(min_stance_s * fs), 1)
    min_gap_n    = max(int(min_gap_s * fs), 1)

    # drop ON (stance) runs shorter than min_stance_n
    lbl, n = label(on)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < min_stance_n:
            on[idx] = False
    # fill OFF (swing) gaps shorter than min_gap_n
    lbl, n = label(~on)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < min_gap_n:
            on[idx] = True
    return on


def tukey_envelope(mask, fs, taper_ms=5.0):
    """Per-stance envelope: 1 inside stance, cosine ramp 0->1 / 1->0 over
    `taper_ms` at each stance boundary, 0 in swing."""
    n_taper = max(int(taper_ms / 1000.0 * fs), 1)
    env = mask.astype(float)

    diff = np.diff(mask.astype(int))
    starts = np.where(diff == 1)[0] + 1
    ends   = np.where(diff == -1)[0] + 1
    if mask[0]:
        starts = np.concatenate(([0], starts))
    if mask[-1]:
        ends = np.concatenate((ends, [len(mask)]))

    ramp_up = np.linspace(0.0, 1.0, n_taper)
    ramp_dn = np.linspace(1.0, 0.0, n_taper)
    for s, e in zip(starts, ends):
        up = slice(s, min(s + n_taper, e))
        env[up] = np.minimum(env[up], ramp_up[:up.stop - up.start])
        dn = slice(max(e - n_taper, s), e)
        env[dn] = np.minimum(env[dn], ramp_dn[-(dn.stop - dn.start):])
    return env


def stance_durations(mask):
    """Lengths (in samples) of the True runs of a boolean mask."""
    d = np.diff(np.concatenate(([0], mask.astype(int), [0])))
    starts = np.where(d == 1)[0]
    ends   = np.where(d == -1)[0]
    return ends - starts


# ---------------------------------------------------------------------------
# Phase 1: GRF 6 Hz filtering
# ---------------------------------------------------------------------------
def filter_grf():
    print("\n[1/3] Filtering GRF at {:.1f} Hz (zero-phase, 4th-order Butterworth) ..."
          .format(GRF_FILTER_FC))
    header, col_labels, data = read_table(GRF_RAW_MOT)
    n_rows, n_cols = data.shape

    t = data[:, 0]
    fs = 1.0 / np.median(np.diff(t))
    nyq = fs / 2.0
    if GRF_FILTER_FC >= nyq:
        raise ValueError("cutoff {:.1f} Hz >= Nyquist {:.1f} Hz".format(GRF_FILTER_FC, nyq))

    b, a = butter(4, GRF_FILTER_FC / nyq, btype="low")
    filt = np.empty_like(data)
    filt[:, 0] = t                      # time column is NOT filtered
    for c in range(1, n_cols):
        filt[:, c] = filtfilt(b, a, data[:, c])

    # quick sanity statistics inside the analysis window
    mask = (t >= TIME_RANGE[0]) & (t <= TIME_RANGE[1])
    print("   rows: raw={} -> filt={} (unchanged), fs={:.0f} Hz".format(n_rows, filt.shape[0], fs))
    print("   window [{} {}]: max vertical GRF raw={:.1f} N -> filt={:.1f} N".format(
        TIME_RANGE[0], TIME_RANGE[1],
        data[mask, 2].max(), filt[mask, 2].max()))
    removed = filt[mask, 1:] - data[mask, 1:]
    print("   RMS of (filtered - raw) in window = {:.3f} (high-freq content removed)".format(
        np.sqrt(np.mean(removed ** 2))))

    write_table(GRF_FILT_MOT, header, col_labels, filt)

    # build a copy of the external loads XML that points to the filtered file
    with open(GRF_RAW_XML, "r", encoding="utf-8") as f:
        xml = f.read()
    datafile_abs = GRF_FILT_MOT.replace(os.sep, "/")   # forward slashes are safest
    xml2 = re.sub(
        r"(<datafile>\s*).*?(\s*</datafile>)",
        lambda m: m.group(1) + datafile_abs + m.group(2),
        xml,
    )
    with open(GRF_FILT_XML, "w", encoding="utf-8") as f:
        f.write(xml2)

    print("   wrote: {}".format(GRF_FILT_MOT))
    print("   wrote: {}".format(GRF_FILT_XML))

    # overlay plots: raw vs filtered
    plot_grf_comparison(t, data, filt,
                        os.path.join(DATA_DIR, "level_walking_grf_filter_comparison_full.png"),
                        downsample=10)
    plot_grf_comparison(t, data, filt,
                        os.path.join(DATA_DIR, "level_walking_grf_filter_comparison_window.png"),
                        tlim=TIME_RANGE)
    plot_grf_comparison(t, data, filt,
                        os.path.join(DATA_DIR, "level_walking_grf_filter_comparison_detail.png"),
                        tlim=(TIME_RANGE[0], TIME_RANGE[0] + 3.0))


# ---------------------------------------------------------------------------
# Phase 1b: Masked GRF filtering (physical constraints)
# ---------------------------------------------------------------------------
def filter_grf_masked():
    """Masked 6 Hz filtering with physical constraints.

    Per belt (independently):
      1. stance = raw vy > MASKED_THR, debounced (short segments removed,
         short gaps filled);
      2. whole-signal 6 Hz Butterworth filtfilt (same as filter_grf);
      3. force/torque (vx vy vz tx ty tz) x per-stance Tukey envelope
         (MASKED_TAPER_MS cosine ramp at each end);
      4. vy = max(vy, 0); wherever vy==0, vx/vz/tx/ty/tz are set to 0;
      5. COP (px py pz): filtered in stance, raw constant value in swing.

    Writes level_walking_grf_masked.mot + .xml and raw-vs-masked PNGs.
    """
    print("\n[1/3] Masked GRF filtering @ {:.1f} Hz "
          "(thr={} N, taper={} ms total) ..."
          .format(GRF_FILTER_FC, MASKED_THR, 2 * MASKED_TAPER_MS))
    header, col_labels, data = read_table(GRF_RAW_MOT)
    n_rows, n_cols = data.shape

    t = data[:, 0]
    fs = 1.0 / np.median(np.diff(t))
    nyq = fs / 2.0
    if GRF_FILTER_FC >= nyq:
        raise ValueError("cutoff {:.1f} Hz >= Nyquist {:.1f} Hz".format(GRF_FILTER_FC, nyq))
    b, a = butter(4, GRF_FILTER_FC / nyq, btype="low")

    filt = np.empty_like(data)
    filt[:, 0] = t
    for c in range(1, n_cols):
        filt[:, c] = filtfilt(b, a, data[:, c])

    # per-belt layout (0-based cols in `data`): belt1=1..9, belt2=10..18
    # offsets within a belt: vx0 vy1 vz2 px3 py4 pz5 tx6 ty7 tz8
    FORCE_TORQUE = (0, 1, 2, 6, 7, 8)   # vx vy vz tx ty tz
    COP          = (3, 4, 5)            # px py pz
    belts_cfg = ((1, 2), (10, 11))      # (first_col, vy_col)

    masks = []
    for start, vy_col in belts_cfg:
        cols = list(range(start, start + 9))
        mask = debounce_mask(data[:, vy_col] > MASKED_THR, fs,
                             MASKED_MIN_STANCE_S, MASKED_MIN_GAP_S)
        env  = tukey_envelope(mask, fs, MASKED_TAPER_MS)

        for off in FORCE_TORQUE:
            filt[:, cols[off]] *= env
        filt[:, vy_col] = np.maximum(filt[:, vy_col], 0.0)

        zero = filt[:, vy_col] == 0.0
        for off in FORCE_TORQUE:
            if off != 1:                       # vx, vz, tx, ty, tz
                filt[zero, cols[off]] = 0.0
        for off in COP:                        # restore raw COP during swing
            filt[~mask, cols[off]] = data[~mask, cols[off]]

        masks.append((mask, zero))

    # --- sanity statistics inside the analysis window ---
    mwin = (t >= TIME_RANGE[0]) & (t <= TIME_RANGE[1])
    print("   rows: {}, fs={:.0f} Hz".format(n_rows, fs))
    for bi, (start, vy_col) in enumerate(belts_cfg):
        mask, zero = masks[bi]
        cols = list(range(start, start + 9))
        durs = stance_durations(mask) / fs
        swing = (~mask) & mwin
        leak = 0.0
        if swing.any():
            leak = max(float(np.abs(filt[swing, cols[off]]).max())
                       for off in FORCE_TORQUE if off != 1)
        print("   belt{}: {} stance events, dur {:.2f}-{:.2f} s, "
              "min(vy)={:.3f} N, swing max|vx,vz,tx,ty,tz|={:.3g}"
              .format(bi + 1, len(durs),
                      durs.min() if len(durs) else float("nan"),
                      durs.max() if len(durs) else float("nan"),
                      filt[mwin, vy_col].min(), leak))

    write_table(GRF_MASKED_MOT, header, col_labels, filt)

    # external-loads XML pointing to the masked file
    with open(GRF_RAW_XML, "r", encoding="utf-8") as f:
        xml = f.read()
    datafile_abs = GRF_MASKED_MOT.replace(os.sep, "/")
    xml2 = re.sub(
        r"(<datafile>\s*).*?(\s*</datafile>)",
        lambda m: m.group(1) + datafile_abs + m.group(2),
        xml,
    )
    with open(GRF_MASKED_XML, "w", encoding="utf-8") as f:
        f.write(xml2)

    print("   wrote: {}".format(GRF_MASKED_MOT))
    print("   wrote: {}".format(GRF_MASKED_XML))

    # overlay plots: raw vs masked
    plot_grf_comparison(t, data, filt,
                        os.path.join(DATA_DIR, "level_walking_grf_masked_comparison_window.png"),
                        tlim=TIME_RANGE, filt_label="belt2 masked",
                        filt_desc="{} Hz masked".format(int(GRF_FILTER_FC)))
    plot_grf_comparison(t, data, filt,
                        os.path.join(DATA_DIR, "level_walking_grf_masked_comparison_detail.png"),
                        tlim=(TIME_RANGE[0], TIME_RANGE[0] + 3.0),
                        filt_label="belt2 masked",
                        filt_desc="{} Hz masked".format(int(GRF_FILTER_FC)))


# ---------------------------------------------------------------------------
# Phase 2: Inverse Kinematics
# ---------------------------------------------------------------------------
def run_ik():
    print("\n[2/3] Running Inverse Kinematics ...")
    tool = osim.InverseKinematicsTool(IK_SETUP)
    tool.set_model_file(MODEL)
    tool.set_marker_file(TRC)
    tool.set_output_motion_file(IK_OUT)
    tool.setStartTime(TIME_RANGE[0])
    tool.setEndTime(TIME_RANGE[1])
    tool.run()

    if not os.path.exists(IK_OUT):
        raise RuntimeError("IK did not produce output: " + IK_OUT)
    _, _, data = read_table(IK_OUT)
    print("   IK output: {} rows, time {:.2f} - {:.2f} s".format(data.shape[0], data[0, 0], data[-1, 0]))
    errs = glob.glob(os.path.join(DATA_DIR, "*marker_errors*.sto"))
    if errs:
        print("   marker error report: {}".format(errs[0]))
    else:
        print("   (no marker error report file found)")


# ---------------------------------------------------------------------------
# Phase 3: Inverse Dynamics (uses filtered GRF)
# ---------------------------------------------------------------------------
def run_id():
    print("\n[3/3] Running Inverse Dynamics (filtered GRF) ...")
    os.makedirs(ID_RESULTS, exist_ok=True)
    tool = osim.InverseDynamicsTool(ID_SETUP)
    tool.setModelFileName(MODEL)
    tool.setCoordinatesFileName(IK_OUT)
    tool.setExternalLoadsFileName(GRF_FILT_XML)
    tool.setResultsDir(ID_RESULTS)
    tool.setStartTime(TIME_RANGE[0])
    tool.setEndTime(TIME_RANGE[1])
    tool.run()

    if not os.path.exists(ID_OUT):
        raise RuntimeError("ID did not produce output: " + ID_OUT)
    header, col_labels, data = read_table(ID_OUT)
    cols = column_names(col_labels)
    if len(cols) != data.shape[1]:
        cols = ["time"] + ["col{}".format(i + 1) for i in range(data.shape[1] - 1)]
    idx = {name: i for i, name in enumerate(cols)}

    mask = (data[:, 0] >= TIME_RANGE[0]) & (data[:, 0] <= TIME_RANGE[1])

    def maxabs_in_window(name):
        return np.max(np.abs(data[mask, idx[name]]))

    print("   ID output: {} rows x {} cols".format(data.shape[0], data.shape[1]))
    print("   [sanity] window [{:.0f} {:.0f}] s:".format(*TIME_RANGE))
    print("     max|pelvis_tx_force| = {:.1f} N".format(maxabs_in_window("pelvis_tx_force")))
    print("     max|pelvis_ty_force| = {:.1f} N".format(maxabs_in_window("pelvis_ty_force")))
    for j in ["hip_flexion_r_moment", "hip_flexion_l_moment",
              "knee_extension_r_moment", "knee_extension_l_moment",
              "ankle_dorsiflexion_r_moment", "ankle_dorsiflexion_l_moment"]:
        if j in idx:
            print("     max|{}| = {:.1f} Nm".format(j, maxabs_in_window(j)))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="OpenSim IK + ID with 6 Hz GRF pre-filtering")
    parser.add_argument("--skip-filter", action="store_true", help="skip GRF filtering")
    parser.add_argument("--skip-ik",    action="store_true", help="skip inverse kinematics")
    parser.add_argument("--skip-id",    action="store_true", help="skip inverse dynamics")
    parser.add_argument("--masked",     action="store_true",
                        help="masked GRF filtering: stance-mask + Tukey ramp + physical constraints")
    args = parser.parse_args()

    global GRF_FILT_MOT, GRF_FILT_XML
    if args.masked:
        GRF_FILT_MOT = GRF_MASKED_MOT
        GRF_FILT_XML = GRF_MASKED_XML

    print("OpenSim version :", osim.GetVersion())
    print("Script          :", os.path.abspath(__file__))
    print("Data dir        :", DATA_DIR)
    check_inputs()

    # run from the data directory so OpenSim resolves any relative paths as in the GUI
    os.chdir(DATA_DIR)

    if not args.skip_filter:
        if args.masked:
            filter_grf_masked()
        else:
            filter_grf()
    if not args.skip_ik:
        backup_if_exists(IK_OUT)
        run_ik()
    if not args.skip_id:
        backup_if_exists(ID_OUT)
        run_id()

    print("\nPipeline finished.")


if __name__ == "__main__":
    main()
