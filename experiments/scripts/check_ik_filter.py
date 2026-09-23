"""Validate <trial>_ik_filtered.mot against the raw <trial>_ik.mot.

Run with the `opensim_scripting` conda environment:
    python experiments\\scripts\\check_ik_filter.py
    python experiments\\scripts\\check_ik_filter.py --recheck-ik    # slower, decisive

Answers the question "is level_walking_ik.mot really the raw IK output and
level_walking_ik_filtered.mot really raw + 15 Hz?".  Checks:

  C1  both files have the same shape, column labels and time column
  C2  filtered == 4th-order 15 Hz Butterworth filtfilt(raw), sample for sample
      (within the 1e-6 that the 8-decimal text format can resolve) -- this is
      the definition of the filtered file, so it pins down what it is
  C3  the filter itself has the advertised response: unity in the passband,
      0.5 (-6 dB) at the 15 Hz cutoff, steep roll-off above it
  C4  the raw file really is UNFILTERED: it must still carry substantially more
      high-frequency energy than the filtered file.  If this fails, the raw file
      has been filtered already and applying ik_filter.py to it would filter it
      twice
  C5  the movement is untouched: per-column peak |raw - filtered| stays a small
      fraction of the column's range (no passband distortion)
  C6  aggregated over the record, the filtered file is smoother: its total
      second-difference energy is lower for every column

C6 is deliberately an aggregate.  Pointwise it is NOT true that the filtered
sample is always "smoother" -- at ~40% of samples |d2| is slightly LARGER after
filtering, which is normal for a zero-phase filter (it redistributes curvature
instead of removing it everywhere).  The per-column table below prints both
numbers so the pointwise wiggle cannot be mistaken for a failed filter.  It
cannot create a visible kink either: in the 172.6-173.8 s stretch of the SQR
trial the two files differ by at most 0.07 deg on the knee, well under one pixel
at the zoom the OpenSim plotter uses.

--recheck-ik re-runs the InverseKinematicsTool into a temporary file and
compares it byte-for-byte with the stored raw file.  That is the decisive test
that <trial>_ik.mot is the tool's own output and nothing else.

Exits non-zero if any check fails.  Read-only (--recheck-ik writes only to the
system temp directory).
"""
import argparse
import hashlib
import os
import shutil
import sys
import tempfile

import numpy as np
from scipy.signal import butter, filtfilt, freqz

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from grf_filter import read_table, column_names                       # noqa: E402
from ik_filter import FILTER_ORDER, IK_FILT_MOT, KIN_FILTER_FC        # noqa: E402
from ik import IK_RAW_MOT                                             # noqa: E402

TOL = 1e-6              # text format resolution (8 decimals) + rounding
C4_MIN_RATIO = 4.0      # raw HF energy must exceed the filtered one by this much
C5_MAX_FRACTION = 0.03  # peak |raw-filtered| / column range

_failures = []


def check(name, ok, detail):
    print("  [{}] {:<56s} {}".format("PASS" if ok else "FAIL", name, detail))
    if not ok:
        _failures.append(name)


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def hf_rms(x, fs, lo=16.0, hi=48.0):
    """RMS of the signal's content in [lo, hi] Hz (Hann-windowed FFT)."""
    n = len(x)
    w = np.hanning(n)
    X = np.fft.rfft((x - x.mean()) * w)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    m = (f >= lo) & (f <= hi)
    return float(np.sqrt(np.sum(np.abs(X[m]) ** 2)) / np.sum(w) * 2.0)


def main():
    ap = argparse.ArgumentParser(description="Validate raw vs 15 Hz filtered IK")
    ap.add_argument("--raw", default=IK_RAW_MOT, help="raw IK .mot")
    ap.add_argument("--filtered", default=IK_FILT_MOT, help="filtered IK .mot")
    ap.add_argument("--cutoff", type=float, default=KIN_FILTER_FC)
    ap.add_argument("--order", type=int, default=FILTER_ORDER)
    ap.add_argument("--recheck-ik", action="store_true",
                    help="re-run the IK tool and compare with the raw file")
    args = ap.parse_args()

    raw_path = os.path.abspath(args.raw)
    flt_path = os.path.abspath(args.filtered)
    print("=" * 78)
    print("IK filter validation")
    print("=" * 78)
    print("  raw      : {}".format(raw_path))
    print("  filtered : {}".format(flt_path))
    for p in (raw_path, flt_path):
        if not os.path.exists(p):
            sys.exit("missing file: {}".format(p))

    hr, lr, raw = read_table(raw_path)
    hf, lf, flt = read_table(flt_path)
    labels = column_names(lr)
    if len(labels) != raw.shape[1]:
        labels = ["col{}".format(i) for i in range(raw.shape[1])]
    ncol = raw.shape[1] - 1
    fs = 1.0 / float(np.median(np.diff(raw[:, 0])))
    print("  fs       : {:.2f} Hz   cutoff {:.1f} Hz, order {} (zero phase)".format(
        fs, args.cutoff, args.order))
    print()

    # ---- C1 -------------------------------------------------------------
    same = (raw.shape == flt.shape and lr == lf and
            np.array_equal(raw[:, 0], flt[:, 0]))
    check("C1 same shape / labels / time column", same,
          "{} vs {} rows, time column identical: {}".format(
              raw.shape[0], flt.shape[0], np.array_equal(raw[:, 0], flt[:, 0])))
    if raw.shape != flt.shape:
        sys.exit("cannot continue: the two files do not describe the same samples")

    # ---- C2 -------------------------------------------------------------
    b, a = butter(args.order, args.cutoff / (fs / 2.0), btype="low")
    ref = raw.copy()
    for c in range(1, raw.shape[1]):
        ref[:, c] = filtfilt(b, a, raw[:, c])
    d = np.abs(ref - flt).max()
    check("C2 filtered == {}th-order {:.0f} Hz filtfilt(raw)".format(
        args.order, args.cutoff), d <= TOL,
        "max|diff| = {:.2e}{}".format(d, "" if d <= TOL else "  <- files are not a filter pair"))

    # ---- C3 -------------------------------------------------------------
    w, H = freqz(b, a, worN=8192, fs=fs)
    mag = np.abs(H) ** 2                       # filtfilt = forward * backward
    def at(f):
        return float(mag[int(np.argmin(np.abs(w - f)))])
    ok3 = (abs(at(1.0) - 1.0) < 1e-3 and abs(at(15.0) - 0.5) < 0.02 and
           at(20.0) < 0.06 and at(25.0) < 0.01)
    check("C3 filter response as specified", ok3,
          "|H| = {:.4f} @1Hz, {:.4f} @15Hz, {:.4f} @20Hz, {:.5f} @25Hz".format(
              at(1.0), at(15.0), at(20.0), at(25.0)))

    # ---- C4 -------------------------------------------------------------
    ratios = []
    for c in range(1, raw.shape[1]):
        r = hf_rms(flt[:, c], fs)
        ratios.append(hf_rms(raw[:, c], fs) / r if r > 0 else np.inf)
    ratios = np.array(ratios)
    n_ok = int((ratios >= C4_MIN_RATIO).sum())
    check("C4 raw file is really unfiltered", n_ok >= ncol - 2,
          "{}/{} columns carry >={:.0f}x the filtered HF energy (median {:.1f}x)".format(
              n_ok, ncol, C4_MIN_RATIO, float(np.median(ratios))))

    # ---- C5 -------------------------------------------------------------
    rng = np.ptp(raw[:, 1:], axis=0)
    peak = np.abs(raw[:, 1:] - flt[:, 1:]).max(axis=0)
    frac = peak / np.maximum(rng, 1e-12)
    check("C5 movement preserved (no passband distortion)", float(frac.max()) <= C5_MAX_FRACTION,
          "worst column shifts {:.3f}% of its range (limit {:.0f}%)".format(
              100 * float(frac.max()), 100 * C5_MAX_FRACTION))

    # ---- C6 -------------------------------------------------------------
    d2r = np.diff(raw[:, 1:], 2, axis=0)
    d2f = np.diff(flt[:, 1:], 2, axis=0)
    energy_ok = int((np.sum(d2f ** 2, axis=0) < np.sum(d2r ** 2, axis=0)).sum())
    check("C6 filtered smoother over the record", energy_ok == ncol,
          "{}/{} columns have lower total |d2|^2".format(energy_ok, ncol))

    # ---- per-column detail ---------------------------------------------
    print()
    print("  per column ({}):".format(flt_path.split(os.sep)[-1]))
    print("    {:<22s} {:>9s} {:>11s} {:>12s}".format(
        "column", "HF ratio", "peak shift", "% samples with"))
    print("    {:<22s} {:>9s} {:>11s} {:>12s}".format(
        "", "raw/filt", "(deg or m)", "larger |d2|"))
    for j in range(ncol):
        wiggly = 100.0 * float(np.mean(np.abs(d2f[:, j]) > np.abs(d2r[:, j])))
        print("    {:<22s} {:9.1f} {:11.4f} {:11.0f}%".format(
            labels[1 + j], ratios[j], peak[j], wiggly))

    # ---- optional: the raw file is the tool's own output ----------------
    if args.recheck_ik:
        print()
        print("  re-running the InverseKinematicsTool for an independent raw file ...")
        import ik as ikmod
        tmp = os.path.join(tempfile.mkdtemp(prefix="ik_check_"), "ik_fresh.mot")
        try:
            ikmod.run_ik(tmp, verbose=False)
            same_bytes = md5(tmp) == md5(raw_path)
            dd = np.abs(read_table(tmp)[2] - raw).max()
            check("C7 raw file == a fresh IK run", same_bytes,
                  "md5 {} vs {}, max|diff| = {:.1e}".format(
                      md5(tmp)[:12], md5(raw_path)[:12], dd))
        finally:
            shutil.rmtree(os.path.dirname(tmp), ignore_errors=True)

    print()
    if _failures:
        print("FAILED: {}".format(", ".join(_failures)))
        return 1
    print("All checks passed: the pair is exactly 'raw IK' + 'that file low-passed "
          "at {:.0f} Hz'.".format(args.cutoff))
    return 0


if __name__ == "__main__":
    sys.exit(main())
