"""Standalone IK filter: raw IK result -> analysis-ready 15 Hz kinematics.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\ik_filter.py"

    level_walking_ik.mot  ->  level_walking_ik_filtered.mot      (defaults)

This step is deliberately separate from IK and from IK/ID/RRA, exactly as
grf_filter.py is separate from the external loads, so that the filter lives in
one place and every consumer sees the same kinematics.  run_ik_id.py and
run_ik_rra.py write the *raw* IK result and then call filter_ik() from here; both
use the same file names, so there is just one pair of files:

    level_walking_ik.mot  ->  level_walking_ik_filtered.mot

Naming -- the two files differ by exactly this one filter pass
-------------------------------------------------------------
    <trial>_ik.mot            InverseKinematicsTool output, NOT filtered
    <trial>_ik_filtered.mot   that file after this filter

Inverse dynamics, RRA and MocoInverse all consume the _filtered file.  Never
point them at the raw one: they differentiate the coordinates twice, so raw IK
noise turns into enormous joint accelerations.

The filter (formerly inlined in run_ik_id.py / run_ik_rra.py as
`filter_kinematics()`, unchanged since):

    fs   = 1 / median(diff(t))                  sampling rate, from the time column
    b, a = butter(4, 15 / (fs/2), btype="low")  digital Butterworth, cutoff 15 Hz
    out[:, c] = filtfilt(b, a, data[:, c])      every column except time

Details that matter
-------------------
* Zero phase.  filtfilt runs the filter forwards and backwards, so no group delay
  is introduced and the sample times stay valid.  A causal lfilter would shift
  the kinematics by ~(order/2)/fc and misalign it against the GRF.
* Pad length.  filtfilt's default padtype='odd' needs at least
  3*max(len(a), len(b)) = 15 samples; a 70 s trial at 100 Hz has 7001, so the
  full record can be filtered in one pass without chunking.
* Cutoff vs Nyquist.  Coordinates are recorded at 100 Hz, so 15 Hz is well below
  the 50 Hz Nyquist limit.  If a file were sampled at <= 30 Hz the filter is not
  defined; the script then stops instead of silently returning wrong data.
* Units.  The Butterworth passband gain is 1 and the filter is linear, so
  filtering degrees (inDegrees=yes) and filtering radians give the same curve --
  no unit conversion is needed or performed.
* The cutoff must stay equal to grf_filter.GRF_FILTER_FC and to the coordinate
  filter in Setup_InverseDynamics.xml.  Inverse dynamics differentiates the
  coordinates twice and adds them to the external loads; if the two carry
  different high-frequency content the residual does not cancel and the joint
  moments are not meaningful.  All three are 15 Hz.

The input file is never modified.

Usage:
    ik_filter.py                                       # the default file above
    ik_filter.py --input <raw_ik.mot>                  # -> <in>_filtered.mot
    ik_filter.py --input <in.mot> --output <out.mot>
    ik_filter.py --input <in.mot> --cutoff 6 --order 4
    ik_filter.py --quiet                               # no per-column summary
"""
import argparse
import os
import sys

import numpy as np
from scipy.signal import butter, filtfilt

# make grf_filter importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import read_table, write_table, column_names   # noqa: E402
from grf_filter import backup_if_exists                        # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

IK_RAW_MOT = os.path.join(DATA_DIR, "level_walking_ik.mot")     # raw IK, do NOT feed to ID
IK_FILT_MOT = os.path.join(DATA_DIR, "level_walking_ik_filtered.mot")

KIN_FILTER_FC = 15.0    # Hz, must equal grf_filter.GRF_FILTER_FC
FILTER_ORDER = 4        # Butterworth order, zero-phase (filtfilt)


def filtered_path(src):
    """`<dir>/<stem>.mot` -> `<dir>/<stem>_filtered.mot`."""
    stem, ext = os.path.splitext(src)
    return stem + "_filtered" + (ext or ".mot")


# ---------------------------------------------------------------------------
# Filter
# ---------------------------------------------------------------------------
def filter_ik(src=IK_RAW_MOT, dst=None, cutoff=KIN_FILTER_FC, order=FILTER_ORDER,
              backup=False, quiet=False):
    """Low-pass every non-time column of `src` and write the result to `dst`.

    `dst` defaults to filtered_path(src).  Called both by this file's main() and
    by run_ik_id.py / run_ik_rra.py right after they write a raw IK result.
    `backup=True` copies an existing `dst` to `dst.bak` first (the pipelines do
    that; the standalone CLI does not).

    Returns (fs, n_columns_filtered, peak_accel_before, peak_accel_after, units).
    """
    if dst is None:
        dst = filtered_path(src)
    if not os.path.exists(src):
        sys.exit("input file not found: {}".format(src))
    if os.path.abspath(dst) == os.path.abspath(src):
        sys.exit("refusing to overwrite the input file in place: {}".format(src))

    if backup:
        bak = backup_if_exists(dst)
        if bak:
            print("   backed up existing file -> {}".format(os.path.basename(bak)))

    header, col_labels, data = read_table(src)
    names = column_names(col_labels)
    t = data[:, 0]

    if data.ndim != 2 or data.shape[1] < 2:
        sys.exit("{} has no coordinate columns".format(src))
    if not np.all(np.isfinite(data)):
        sys.exit("{} contains NaN/Inf -- fix the IK result first".format(src))
    if np.any(np.diff(t) <= 0):
        sys.exit("{} has a non-increasing time column".format(src))

    fs = 1.0 / float(np.median(np.diff(t)))
    nyq = fs / 2.0
    print("   {} -> {}".format(os.path.basename(src), os.path.basename(dst)))
    print("   input : {} rows x {} cols, {:.2f}-{:.2f} s".format(
        data.shape[0], data.shape[1], t[0], t[-1]))
    print("   fs    : {:.2f} Hz  (dt = {:.6f} s, median of {} steps)".format(
        fs, 1.0 / fs, data.shape[0] - 1))

    if cutoff >= nyq:
        sys.exit("cutoff {:.1f} Hz >= Nyquist {:.1f} Hz (fs = {:.1f} Hz): "
                 "the filter is not defined for this file".format(cutoff, nyq, fs))

    b, a = butter(order, cutoff / nyq, btype="low")
    padlen = 3 * max(len(a), len(b))
    if data.shape[0] <= padlen:
        sys.exit("{} rows is too few for a {}th-order filtfilt (needs > {} samples)"
                 .format(data.shape[0], order, padlen))

    out = data.copy()
    for c in range(1, data.shape[1]):
        out[:, c] = filtfilt(b, a, data[:, c])

    write_table(dst, header, col_labels, out)

    # a filtered coordinate must keep its mean and lose its high-frequency spikes.
    # axis=0 is essential: without it np.gradient also differences across columns,
    # i.e. between unrelated coordinates, and reports a meaningless number.
    dt = 1.0 / fs
    acc_before = np.abs(np.gradient(np.gradient(data[:, 1:], dt, axis=0), dt, axis=0))
    acc_after = np.abs(np.gradient(np.gradient(out[:, 1:], dt, axis=0), dt, axis=0))
    units = "deg" if any("inDegrees=yes" in h for h in header) else "rad"

    print("   filtered {} coordinate columns at {:.1f} Hz, order {} "
          "(fs = {:.0f} Hz)".format(data.shape[1] - 1, cutoff, order, fs))
    print("   peak |accel| over all columns: {:.0f} -> {:.0f} {}/s^2".format(
        float(acc_before.max()), float(acc_after.max()), units))
    report(data, out, names, quiet=quiet)
    print("   wrote: {}".format(dst))

    return (fs, data.shape[1] - 1, float(acc_before.max()),
            float(acc_after.max()), units)


def report(raw, flt, names, quiet=False):
    """Quantify what the filter removed, per column (printed unless quiet)."""
    if len(names) != raw.shape[1]:
        names = ["col{}".format(i) for i in range(raw.shape[1])]

    band = raw[:, 1:] - flt[:, 1:]                     # what was removed
    peak_raw = np.abs(raw[:, 1:]).max(axis=0)
    peak_rm = np.abs(band).max(axis=0)
    rms_raw = raw[:, 1:].std(axis=0)
    rms_rm = band.std(axis=0)
    dc = np.abs(flt[:, 1:].mean(axis=0) - raw[:, 1:].mean(axis=0))
    frac = 100.0 * rms_rm / np.maximum(rms_raw, 1e-12)

    if not quiet:
        print("   per-column summary (units: deg for rotations, m for translations):")
        print("      {:<22s} {:>10s} {:>10s} {:>8s} {:>10s}".format(
            "column", "peak", "removed", "rem%", "DC shift"))
        for j, nm in enumerate(names[1:]):
            print("      {:<22s} {:10.4f} {:10.5f} {:7.2f}% {:10.2e}".format(
                nm, peak_raw[j], peak_rm[j], frac[j], dc[j]))

    j = int(np.argmax(frac))
    print("   removed {:.2f}% of the signal RMS (median over columns), "
          "largest: {} ({:.2f}%)".format(
              float(np.median(frac)), names[1 + j], frac[j]))
    print("   total DC shift over all columns: {:.2e} (expect ~0: filtfilt is "
          "zero-phase and unity-gain at DC)".format(float(dc.max())))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Raw IK .mot -> {} Hz, order {} Butterworth (zero phase)"
                    .format(KIN_FILTER_FC, FILTER_ORDER))
    ap.add_argument("--input", default=IK_RAW_MOT,
                    help="raw (unfiltered) IK .mot to read")
    ap.add_argument("--output", default=None,
                    help="filtered .mot to write (default: <input>_filtered.mot)")
    ap.add_argument("--cutoff", type=float, default=KIN_FILTER_FC,
                    help="cutoff in Hz (default {})".format(KIN_FILTER_FC))
    ap.add_argument("--order", type=int, default=FILTER_ORDER,
                    help="Butterworth order (default {})".format(FILTER_ORDER))
    ap.add_argument("--quiet", action="store_true", help="less per-column output")
    ap.add_argument("--backup", action="store_true",
                    help="copy an existing output to <output>.bak first")
    args = ap.parse_args()

    print("IK filter: raw -> {:.1f} Hz, order {}, zero-phase (filtfilt)".format(
        args.cutoff, args.order))
    filter_ik(os.path.abspath(args.input), args.output, args.cutoff, args.order,
              backup=args.backup, quiet=args.quiet)
    print("\nDone.  The input file was not modified.")


if __name__ == "__main__":
    main()
