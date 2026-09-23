"""Full-time GRF filter: raw record -> 6 Hz low-passed external loads.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\grf_filter_fulltime.py"

    level_walking_grf.mot  ->  level_walking_grf_filtered_6Hz_fulltime.mot
    level_walking_grf.xml  ->  level_walking_grf_filtered_6Hz_fulltime.xml

What it does, and how it differs from grf_filter.py
--------------------------------------------------
grf_filter.py detects the stance phases first and filters each one only, with
mirror padding, a Tukey taper and the swing phase forced to exactly zero; it also
leaves the COP columns untouched.  This script deliberately does none of that:

    * all 9 columns per belt (force vx vy vz, COP px py pz, moment tx ty tz) are
      low-pass filtered, COP included;
    * the whole record is filtered in one pass -- every sample, swing included,
      with no contact detection, no taper and no padding;
    * nothing is forced to zero, so the swing phase is whatever the filter leaves
      there rather than exactly 0.

The recipe is the project's usual one, so the result is band-compatible with the
kinematics and with Setup_InverseDynamics.xml:

    b, a = butter(4, fc / (fs/2), btype="low");  y = filtfilt(b, a, x)

Trade-off to be aware of: a zero-phase filter run over the whole record smears
every step in it over roughly +/-4/fc seconds -- at 6 Hz that is about +/-0.67 s.
The plate's load/unload steps at touch-down and toe-off and the swing-to-stance
COP jump therefore bleed into the first and last ~0.7 s of each stance phase.
The stance-segmented grf_filter.py exists precisely to avoid that; this script
exists to give the plain, whole-record alternative for comparison, and it prints
the size of that bleed (swing-phase force leak, COP step loss) so the cost is
visible rather than assumed.

Files default to the SQR level-walking trial.  The cutoff is a module constant
with a --cutoff override; the output path follows --output or defaults to
<trial>_grf_filtered_<fc>Hz_fulltime.mot, so this script cannot silently
overwrite grf_filter.py's level_walking_grf_filtered.mot.
"""
import argparse
import os
import sys

import numpy as np
from scipy.signal import butter, filtfilt

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import (read_table, write_table, column_names, write_loads_xml,
                        debounce_mask, BELTS, OFF_FORCE, OFF_COP, OFF_TORQUE,
                        FILT_THR)                                  # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

GRF_RAW_MOT = os.path.join(DATA_DIR, "level_walking_grf.mot")
GRF_RAW_XML = os.path.join(DATA_DIR, "level_walking_grf.xml")

GRF_FILTER_FC = 6.0         # Hz -- this script's whole point
FILTER_ORDER = 4            # Butterworth order, zero-phase (filtfilt)

# Reporting only: the analysis window the rest of the project works in.
TIME_RANGE = (170.0, 240.0)

COLUMN_GROUPS = (("force", OFF_FORCE), ("COP", OFF_COP), ("moment", OFF_TORQUE))


def default_output(cutoff=GRF_FILTER_FC):
    """<trial>_grf.mot -> <trial>_grf_filtered_<fc>Hz_fulltime.mot (+ .xml)."""
    tag = "{:g}Hz".format(cutoff)
    stem = os.path.basename(GRF_RAW_MOT).replace("_grf.mot", "_grf_filtered_{}_fulltime.mot".format(tag))
    return os.path.join(DATA_DIR, stem), os.path.join(DATA_DIR, stem[:-4] + ".xml")


# ---------------------------------------------------------------------------
# Filter
# ---------------------------------------------------------------------------
def filter_grf_fulltime(cutoff=GRF_FILTER_FC, order=FILTER_ORDER,
                        out_mot=None, out_xml=None, quiet=False):
    """Low-pass all 18 data columns of the raw GRF over the whole record.

    Returns (out_mot, out_xml, info) where info is a dict of the reported
    numbers, so a caller can assert on them.
    """
    out_mot = out_mot or default_output(cutoff)[0]
    out_xml = out_xml or default_output(cutoff)[1]

    print("GRF full-time filter: raw -> {:.1f} Hz, order {}, whole record, "
          "no stance segmentation".format(cutoff, order))
    header, col_labels, data = read_table(GRF_RAW_MOT)
    names = column_names(col_labels)

    t = data[:, 0]
    fs = 1.0 / float(np.median(np.diff(t)))
    nyq = fs / 2.0
    if cutoff >= nyq:
        sys.exit("cutoff {:.1f} Hz >= Nyquist {:.1f} Hz".format(cutoff, nyq))
    b, a = butter(order, cutoff / nyq, btype="low")
    print("   rows {}, fs {:.0f} Hz, cutoff {:.1f} Hz, window {:.1f}-{:.1f} s".format(
        data.shape[0], fs, cutoff, t[0], t[-1]))

    # verify the assumed column layout against the header, as grf_filter.py does
    if names and len(names) == data.shape[1]:
        for bi, (start, vy_col) in enumerate(BELTS):
            exp_vy = "ground_force{}_vy".format(bi + 1)
            if names[vy_col] != exp_vy:
                sys.exit("column layout mismatch: column {} is '{}', expected '{}'"
                         .format(vy_col, names[vy_col], exp_vy))
        print("   header layout verified against BELTS/OFF_* constants")

    # every non-time column, in one pass over the whole record
    out = data.copy()
    for c in range(1, data.shape[1]):
        out[:, c] = filtfilt(b, a, data[:, c])

    info = {"fs": fs, "cutoff": cutoff, "order": order,
            "removed_pct": {}, "swing_leak": {}, "cop_step": {}}

    print("\n   per column (all 18 filtered; removed RMS as a fraction of the "
          "signal RMS):")
    for bi, (start, vy_col) in enumerate(BELTS):
        print("      belt{}:".format(bi + 1))
        for grp, offs in COLUMN_GROUPS:
            for o in offs:
                c = start + o          # BELTS holds the block's first column
                nm = names[c] if c < len(names) else "col{}".format(c)
                short = nm.replace("ground_force1_", "").replace("ground_force2_", "")
                rem = 100.0 * (data[:, c] - out[:, c]).std() / max(data[:, c].std(), 1e-12)
                info["removed_pct"][nm] = rem
                print("         {:<5s} {:<6s} removed {:5.1f}%   peak {:10.4f} -> {:10.4f}".format(
                    short, grp, rem, float(np.abs(data[:, c]).max()),
                    float(np.abs(out[:, c]).max())))

    # ---- what the whole-record filter costs -----------------------------
    # The plate reads a constant (zero force, plate-centre COP) while unloaded, so
    # any non-zero value in swing is bleed from the neighbouring stance phase.
    mwin = (t >= TIME_RANGE[0]) & (t <= TIME_RANGE[1])
    print("\n   cost of filtering the whole record (6 Hz: +/-{:.2f} s smearing):".format(
        4.0 / cutoff))
    for bi, (start, vy_col) in enumerate(BELTS):
        mask = debounce_mask(data[:, vy_col] > FILT_THR, fs)
        swing = (~mask) & mwin
        stance = mask & mwin
        leak_raw = float(np.abs(data[swing, vy_col]).max()) if swing.any() else 0.0
        leak_fil = float(np.abs(out[swing, vy_col]).max()) if swing.any() else 0.0
        pk_raw = float(data[stance, vy_col].max())
        pk_fil = float(out[stance, vy_col].max())
        info["swing_leak"]["belt{}".format(bi + 1)] = leak_fil

        # swing -> stance COP step (a real measurement; filtering smooths it)
        cop_col = start + OFF_COP[0]
        on = np.where(np.diff(mask.astype(int)) == 1)[0] + 1
        on = on[mwin[on]] if on.size else on
        steps_r = [abs(data[i, cop_col] - data[i - 1, cop_col]) for i in on]
        steps_f = [abs(out[i, cop_col] - out[i - 1, cop_col]) for i in on]
        step_r = float(np.median(steps_r)) if steps_r else float("nan")
        step_f = float(np.median(steps_f)) if steps_f else float("nan")
        info["cop_step"]["belt{}".format(bi + 1)] = (step_r, step_f)

        print("      belt{}: vy peak {:8.2f} -> {:8.2f} N (kept {:.2f}%)".format(
            bi + 1, pk_raw, pk_fil, 100.0 * pk_fil / pk_raw))
        print("             swing |vy| raw {:.3f} -> filtered {:.3f} N "
              "(non-zero = bleed)".format(leak_raw, leak_fil))
        print("             swing->stance COP step median {:.4f} -> {:.4f} m "
              "(kept {:.0f}%)".format(step_r, step_f, 100.0 * step_f / step_r))

    write_table(out_mot, header, col_labels, out)
    write_loads_xml(out_xml, out_mot, template_xml=GRF_RAW_XML)
    print("\n   wrote: {}".format(out_mot))
    print("   wrote: {}".format(out_xml))
    if not quiet:
        print("\n   NOTE: the swing phase is NOT forced to zero here (grf_filter.py "
              "does that).\n         Use level_walking_grf_filtered.xml (15 Hz, "
              "stance-segmented) for\n         inverse dynamics unless you have a "
              "reason to prefer this one.")
    return out_mot, out_xml, info


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Raw GRF -> full-time {:.0f} Hz Butterworth filter (no stance "
                    "segmentation)".format(GRF_FILTER_FC))
    ap.add_argument("--cutoff", type=float, default=GRF_FILTER_FC,
                    help="cutoff in Hz (default {})".format(GRF_FILTER_FC))
    ap.add_argument("--order", type=int, default=FILTER_ORDER,
                    help="Butterworth order (default {})".format(FILTER_ORDER))
    ap.add_argument("--output", default=None,
                    help="output .mot (default: <trial>_grf_filtered_<fc>Hz_fulltime.mot)")
    ap.add_argument("--quiet", action="store_true", help="less output")
    args = ap.parse_args()

    if not os.path.exists(GRF_RAW_MOT):
        sys.exit("missing input: {}".format(GRF_RAW_MOT))
    mot, xml, _ = filter_grf_fulltime(args.cutoff, args.order, out_mot=args.output,
                                      quiet=args.quiet)
    print("\nDone.  External loads XML: {}".format(xml))


if __name__ == "__main__":
    main()
