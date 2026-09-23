"""Standalone GRF filter: raw record -> analysis-ready external loads.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\grf_filter.py"

    level_walking_grf.mot  ->  level_walking_grf_filtered.mot
    level_walking_grf.xml  ->  level_walking_grf_filtered.xml

This step is deliberately separate from IK/ID and from RRA so that both pipelines
consume exactly the same external loads and cannot drift apart.

Per belt, independently:

  1. stance = raw vy > FILT_THR, debounced (short runs dropped, short gaps
     filled);
  2. force and free-moment columns (vx vy vz tx ty tz) are low-pass filtered
     *inside each stance segment only*, with mirror padding, so the filter never
     sees the load/unload step at a segment end.  Swing is exactly zero;
  3. those columns are then multiplied by a Tukey envelope (FILT_TAPER_MS cosine
     ramp at each end of every stance phase).  The ramp sits on the low-force
     ends of the stance phase, so it removes the plate's load/unload steps
     without touching the peak force;
  4. the horizontal torque columns (tx, tz) are set to exactly 0: in this record
     they hold nothing but ~1e-5 Nm rounding noise, while ty -- the free moment
     about the vertical axis -- carries real signal and is filtered normally;
  5. the COP columns (px, py, pz) are filtered the SAME segment-wise way -- inside
     each stance phase, with mirror padding -- and left untouched during swing.

Why the COP is filtered segment-wise and never across the swing/stance boundary
------------------------------------------------------------------------------
While the plate is unloaded it reports its geometric centre as a constant
placeholder, so at touch-down the COP jumps ~0.5 m to the real contact point
within one sample.  A zero-phase filter applied to the whole record reads that
step as signal and rings across about +/-4/fc seconds, corrupting the first
~100 ms of every stance phase (0.10 m RMS, 0.65 m peak in this record) and
injecting up to 31 Nm of spurious moment into the loading response.

Filtering each stance phase on its own, mirror-padded at both ends, avoids that
coupling entirely: the filter never sees the step.  Two COP-specific details:

  * swing is NOT zeroed -- the placeholder is the value OpenSim must read while
    the force is zero, and it is a valid position;
  * the segment ends are left to the mirror padding.  Pinning the filtered
    trajectory back onto the raw first sample was tried and is much worse: that
    sample can be a COP outlier (the plate divides by a near-zero Fy at contact),
    and a constant shift moves the whole stance trajectory -- 0.14 m RMS COP
    error and up to 95 Nm of spurious moment, versus 0.014 m and 1.6 Nm when the
    ends are left alone.  The small remaining swing->stance step is harmless
    because the force is only 20-40 N at that instant.

Filtering the COP is what makes the external-load set dynamically consistent:
the moment OpenSim applies is r x F, so a COP and a force filtered at the same
cutoff belong together.

Usage:
    grf_filter.py            # writes level_walking_grf_filtered.mot / .xml
    grf_filter.py --quiet    # same, less chatter

Verify the result with experiments/scripts/check_grf_cop.py
"""
import argparse
import os
import re
import shutil
import sys

import numpy as np
from scipy.ndimage import label
from scipy.signal import butter, filtfilt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

GRF_RAW_MOT = os.path.join(DATA_DIR, "level_walking_grf.mot")
GRF_RAW_XML = os.path.join(DATA_DIR, "level_walking_grf.xml")
GRF_FILT_MOT = os.path.join(DATA_DIR, "level_walking_grf_filtered.mot")
GRF_FILT_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")

# The single cutoff shared by the GRF, the IK coordinates and
# Setup_InverseDynamics.xml.  If they disagree the external loads and the
# kinematics are band-limited differently and inverse dynamics / MocoInverse are
# not self-consistent.
GRF_FILTER_FC = 15.0    # Hz
FILTER_ORDER = 4        # Butterworth order, zero-phase (filtfilt)
TIME_RANGE = (170.0, 240.0)     # analysis window, for reporting only

# contact detection / debounce / taper
FILT_THR = 10.0             # N, vertical-force contact threshold
FILT_MIN_STANCE_S = 0.10    # drop stance segments shorter than this
FILT_MIN_GAP_S = 0.04       # fill swing gaps shorter than this
FILT_TAPER_MS = 25.0        # per-end cosine ramp (50 ms total)
FILT_PAD_S = 0.30           # mirror padding per stance segment before filtering

# Per-belt column layout inside the data array (column 0 is time):
#   cols 1-9   = belt1 vx vy vz | px py pz | tx ty tz
#   cols 10-18 = belt2 vx vy vz | px py pz | tx ty tz
OFF_FORCE = (0, 1, 2)       # vx, vy (vertical), vz
OFF_COP = (3, 4, 5)         # px, py, pz  -- NOT filtered
OFF_TORQUE = (6, 7, 8)      # tx, ty (free moment about vertical), tz
BELTS = ((1, 2), (10, 11))  # (first column of the belt, vertical-force column)


# ---------------------------------------------------------------------------
# OpenSim text-table I/O -- the shared home for the whole script family
# (ik_filter.py, ik.py, id.py, rra.py and both run_ik_*.py import from here)
# ---------------------------------------------------------------------------
def read_table(path):
    """Read an OpenSim .mot/.sto file -> (header_lines, column_labels, ndarray).

    The column-labels line sits right after 'endheader'.
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

    data = np.loadtxt([ln for ln in lines[i:] if ln.strip()])
    return header, col_labels, data


def write_table(path, header, col_labels, data):
    """Write an OpenSim .mot/.sto file, preserving the original header format."""
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
    """Copy `path` to `<path>.bak` if it exists; return the backup path or None.

    Every stage of the pipeline overwrites its previous output in place (as the
    OpenSim GUI does), so each one keeps the pre-run copy under .bak.  The
    caller decides whether to announce it.
    """
    if not os.path.exists(path):
        return None
    bak = path + ".bak"
    shutil.copy2(path, bak)
    return bak


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------
def debounce_mask(on, fs, min_stance_s=FILT_MIN_STANCE_S, min_gap_s=FILT_MIN_GAP_S):
    """Clean a boolean contact mask: drop stance runs shorter than
    `min_stance_s` and fill swing gaps shorter than `min_gap_s`."""
    on = on.copy()
    min_stance_n = max(int(min_stance_s * fs), 1)
    min_gap_n = max(int(min_gap_s * fs), 1)

    lbl, n = label(on)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < min_stance_n:
            on[idx] = False
    lbl, n = label(~on)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < min_gap_n:
            on[idx] = True
    return on


def tukey_envelope(mask, fs, taper_ms=FILT_TAPER_MS):
    """Per-stance envelope: 1 inside stance, cosine ramp 1->0 / 0->1 over
    `taper_ms` at each stance boundary, 0 in swing."""
    n_taper = max(int(taper_ms / 1000.0 * fs), 1)
    env = mask.astype(float)

    diff = np.diff(mask.astype(int))
    starts = np.where(diff == 1)[0] + 1
    ends = np.where(diff == -1)[0] + 1
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
    ends = np.where(d == -1)[0]
    return ends - starts


def segment_filter(x, mask, b, a, fs, pad_s=FILT_PAD_S):
    """Low-pass `x` one stance segment at a time, with mirror padding.

    The raw record is a concatenation of physically unrelated pieces: while the
    plate is unloaded its COP holds the plate-centre constant and its
    force/torque channels hold sensor noise.  Filtering the whole record in one
    pass therefore lets a zero-phase filter smear that unloaded stretch across
    the beginning of every stance phase -- exactly the loading response.
    Filtering each stance run separately, with its ends mirror-reflected so the
    filter sees no step, removes that coupling.

    Returns an array equal to `x` filtered inside stance and exactly 0 outside.
    """
    out = np.zeros_like(x)
    pad = max(int(round(pad_s * fs)), 1)
    lbl, n = label(mask)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < 3:                      # too short to pad either side
            out[idx] = x[idx]
            continue
        s, e = idx[0], idx[-1] + 1
        seg = x[s:e]
        m = min(pad, len(seg) - 1)
        ext = np.concatenate((seg[1:m + 1][::-1], seg, seg[-m - 1:-1][::-1]))
        y = filtfilt(b, a, ext)
        out[s:e] = y[m:m + len(seg)]
    return out


def cop_segment_filter(x, mask, b, a, fs, pad_s=FILT_PAD_S):
    """Low-pass the COP one stance segment at a time.

    Same mirror-padded segment filtering as `segment_filter`, with one COP
    specific difference: swing is NOT zeroed.  During swing the COP column holds
    the plate's geometric centre -- a placeholder, but a valid position, and the
    value OpenSim must read while the force is zero -- so it is left exactly as
    measured.

    Returns an array that is `x` filtered inside stance and equal to `x` in swing.

    Note on the segment ends: a mirror-padded segment filter pulls the two end
    samples a few mm toward the segment interior, so the filtered value at
    touch-down is not bit-identical to the measured contact point.  That is
    deliberate -- pinning it back was tried and is worse.  The raw first stance
    sample can be a COP outlier (the plate divides by a near-zero Fy at contact,
    so px has been seen to jump 0.65 m for one or two samples); a constant shift
    that forces that sample to match moves the WHOLE stance trajectory, giving a
    0.14 m RMS COP error and up to 95 Nm of spurious moment.  Leaving the ends to
    the mirror padding gives 0.014 m RMS and ~1.6 Nm instead.

    The residual swing->stance step is harmless: r x F with F ~ 20-40 N at
    contact makes the applied moment continuous anyway.
    """
    out = x.copy()
    pad = max(int(round(pad_s * fs)), 1)
    lbl, n = label(mask)
    for i in range(1, n + 1):
        idx = np.where(lbl == i)[0]
        if len(idx) < 3:                      # too short to pad either side
            continue
        s, e = idx[0], idx[-1] + 1
        seg = x[s:e]
        m = min(pad, len(seg) - 1)
        ext = np.concatenate((seg[1:m + 1][::-1], seg, seg[-m - 1:-1][::-1]))
        out[s:e] = filtfilt(b, a, ext)[m:m + len(seg)]
    return out


def write_loads_xml(dst_xml, mot_path, template_xml=GRF_RAW_XML):
    """Copy an ExternalLoads XML and repoint <datafile> at `mot_path`."""
    with open(template_xml, "r", encoding="utf-8") as f:
        xml = f.read()
    datafile_abs = mot_path.replace(os.sep, "/")     # forward slashes are safest
    xml2 = re.sub(r"(<datafile>\s*).*?(\s*</datafile>)",
                  lambda m: m.group(1) + datafile_abs + m.group(2),
                  xml, flags=re.S)
    with open(dst_xml, "w", encoding="utf-8") as f:
        f.write(xml2)


# ---------------------------------------------------------------------------
# Main filter
# ---------------------------------------------------------------------------
def filter_grf(quiet=False):
    """Build level_walking_grf_filtered.mot / .xml from the raw record."""
    print("GRF filter: raw -> {:.1f} Hz, order {}, per-stance mirror-padded "
          "segments".format(GRF_FILTER_FC, FILTER_ORDER))
    header, col_labels, data = read_table(GRF_RAW_MOT)
    names = column_names(col_labels)

    t = data[:, 0]
    fs = 1.0 / np.median(np.diff(t))
    nyq = fs / 2.0
    if GRF_FILTER_FC >= nyq:
        sys.exit("cutoff {:.1f} Hz >= Nyquist {:.1f} Hz".format(GRF_FILTER_FC, nyq))
    b, a = butter(FILTER_ORDER, GRF_FILTER_FC / nyq, btype="low")

    out = data.copy()
    mwin = (t >= TIME_RANGE[0]) & (t <= TIME_RANGE[1])
    print("   rows {}, fs {:.0f} Hz, cutoff {:.1f} Hz".format(
        data.shape[0], fs, GRF_FILTER_FC))

    # verify the assumed column layout against the header
    if names and len(names) == data.shape[1]:
        for bi, (start, vy_col) in enumerate(BELTS):
            exp_vy = "ground_force{}_vy".format(bi + 1)
            if names[vy_col] != exp_vy:
                sys.exit("column layout mismatch: column {} is '{}', expected '{}'"
                         .format(vy_col, names[vy_col], exp_vy))
        print("   header layout verified against BELTS/OFF_* constants")

    for bi, (start, vy_col) in enumerate(BELTS):
        fcols = [start + o for o in OFF_FORCE]
        tcols = [start + o for o in OFF_TORQUE]
        ccols = [start + o for o in OFF_COP]

        mask = debounce_mask(data[:, vy_col] > FILT_THR, fs)
        env = tukey_envelope(mask, fs)

        for c in fcols + tcols:
            out[:, c] = segment_filter(data[:, c], mask, b, a, fs) * env

        out[:, tcols[0]] = 0.0                      # tx: rounding noise only
        out[:, tcols[2]] = 0.0                      # tz: rounding noise only

        # COP: filtered inside stance, untouched in swing -- see
        # cop_segment_filter() for why the ends are left to the padding
        for c in ccols:
            out[:, c] = cop_segment_filter(data[:, c], mask, b, a, fs)

        swing = (~mask) & mwin
        filtered_cols = fcols + tcols
        swing_max = (float(np.abs(out[swing][:, filtered_cols]).max())
                     if swing.any() else 0.0)
        stance = mask & mwin
        raw_pk = float(data[stance, vy_col].max())
        out_pk = float(out[stance, vy_col].max())
        raw_ty = float(np.abs(data[stance, tcols[1]]).max())
        out_ty = float(np.abs(out[stance, tcols[1]]).max())
        durs = stance_durations(mask) / fs

        print("   belt{}: {} stance events, dur {:.2f}-{:.2f} s, duty {:.1f}%".format(
            bi + 1, len(durs),
            durs.min() if len(durs) else float("nan"),
            durs.max() if len(durs) else float("nan"),
            100.0 * mask[mwin].mean()))
        print("          vy peak   raw {:8.2f} -> filtered {:8.2f} N  (kept {:.2f}%)".format(
            raw_pk, out_pk, 100.0 * out_pk / raw_pk))
        print("          |ty| peak raw {:8.3f} -> filtered {:8.3f} Nm (kept {:.1f}%)".format(
            raw_ty, out_ty, 100.0 * out_ty / raw_ty))
        if not quiet:
            print("          swing max|force,torque| = {:.3e} (expect exactly 0)".format(
                swing_max))
            # COP: filtered in stance, identical to raw in swing, and the
            # swing->stance step must remain the measured one
            for c in ccols:
                cn = names[c] if c < len(names) else "col{}".format(c)
                if not np.array_equal(out[swing, c], data[swing, c]):
                    sys.exit("COP column {} was modified in swing".format(cn))
                d_st = float(np.abs(out[stance, c] - data[stance, c]).max())
                print("          COP {:<16s} swing untouched, stance max|diff| "
                      "{:.5f} m".format(cn, d_st))
            on = np.where(np.diff(mask.astype(int)) == 1)[0] + 1
            on = on[mwin[on]] if on.size else on
            on = on[on > 0]
            # the COP error that actually matters is the moment it creates
            dM = np.abs(out[stance][:, ccols] - data[stance][:, ccols]) * \
                data[stance, vy_col][:, None]
            print("          COP moment error in stance: max {:.2f} Nm, "
                  "rms {:.2f} Nm".format(dM.max(), np.sqrt(np.mean(dM ** 2))))
            for c in (ccols[0], ccols[2]):
                if on.size:
                    step_r = np.median([abs(data[o, c] - data[o - 1, c]) for o in on])
                    step_f = np.median([abs(out[o, c] - out[o - 1, c]) for o in on])
                    print("          COP {:<16s} stance-entry step {:.4f} -> "
                          "{:.4f} m (n={})".format(
                              names[c] if c < len(names) else str(c),
                              step_r, step_f, int(on.size)))
        if swing_max != 0.0:
            sys.exit("swing phase is not exactly zero for belt{}".format(bi + 1))

    write_table(GRF_FILT_MOT, header, col_labels, out)
    write_loads_xml(GRF_FILT_XML, GRF_FILT_MOT)
    print("   wrote: {}".format(GRF_FILT_MOT))
    print("   wrote: {}".format(GRF_FILT_XML))


def main():
    ap = argparse.ArgumentParser(description="Raw GRF -> filtered external loads")
    ap.add_argument("--quiet", action="store_true", help="less per-column output")
    args = ap.parse_args()

    if not (os.path.exists(GRF_RAW_MOT) and os.path.exists(GRF_RAW_XML)):
        sys.exit("missing input: {} / {}".format(GRF_RAW_MOT, GRF_RAW_XML))
    filter_grf(quiet=args.quiet)
    print("\nDone.  Validate with: check_grf_cop.py")


if __name__ == "__main__":
    main()
