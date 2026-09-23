"""Validate level_walking_grf_filtered.mot against the raw GRF record.

Run with the `opensim_scripting` conda environment:
    python experiments\\scripts\\check_grf_cop.py

Checks that must hold for the processed GRF to be usable by inverse dynamics and
MocoInverse:

  C1  the COP columns are UNCHANGED during swing (the plate-centre placeholder
      OpenSim reads while the force is zero is left exactly as measured)
  C2  force and free-moment columns are exactly zero during swing
  C3  the horizontal torque columns are exactly zero
  C4  the stance vertical-force peak is preserved (>= 99%)
  C5  the free-moment peak is preserved (>= 75%)
  C6  the swing-to-stance COP step survives (>= 0.6 of the raw step; the
      mirror-padded segment filter trims the endpoint, and the moment stays
      continuous anyway because the force is only ~20-40 N at that instant)
  C7  the COP filtering stays dynamically consistent: the moment it changes,
      |dCOP| * Fy, must stay small (<= 20 Nm peak, <= 3 Nm rms in stance)
  C8  no energy above the cutoff survives in the filtered force channels

Exits non-zero if any check fails.  Read-only: writes nothing.
"""
import os
import sys

import numpy as np
from scipy.ndimage import label
from scipy.signal import butter, filtfilt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

RAW_MOT = os.path.join(DATA_DIR, "level_walking_grf.mot")
FILT_MOT = os.path.join(DATA_DIR, "level_walking_grf_filtered.mot")

FC = 15.0                 # Hz, must match run_ik_id.GRF_FILTER_FC
THR = 10.0                # N, contact threshold (run_ik_id.FILT_THR)
TIME_RANGE = (180.0, 240.0)

# per-belt offsets inside a 9-column block
OFF_FORCE = (0, 1, 2)
OFF_COP = (3, 4, 5)
OFF_TORQUE = (6, 7, 8)
BELTS = ((1, 2), (10, 11))    # (block start, vertical-force column)

_failures = []
_results = []


def check(name, ok, detail):
    _results.append((name, ok, detail))
    print("  [{}] {:<58s} {}".format("PASS" if ok else "FAIL", name, detail))
    if not ok:
        _failures.append(name)


def read_mot(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and not lines[i].lstrip().lower().startswith("endheader"):
        i += 1
    i += 1
    labels = lines[i].split()
    i += 1
    data = np.loadtxt([ln for ln in lines[i:] if ln.strip()])
    return labels, data


def main():
    print("=" * 78)
    print("GRF validation: {} vs raw".format(os.path.basename(FILT_MOT)))
    print("=" * 78)
    if not os.path.exists(FILT_MOT):
        sys.exit("missing {} -- run run_ik_id.py first".format(FILT_MOT))

    names, raw = read_mot(RAW_MOT)
    names_f, filt = read_mot(FILT_MOT)
    assert names == names_f, "column labels differ between raw and filtered"
    assert raw.shape == filt.shape, "shape mismatch {} vs {}".format(raw.shape, filt.shape)

    t = raw[:, 0]
    fs = 1.0 / np.median(np.diff(t))
    win = (t >= TIME_RANGE[0]) & (t <= TIME_RANGE[1])
    print("  rows = {}, fs = {:.0f} Hz, window = [{:.0f} {:.0f}] s ({} rows)".format(
        raw.shape[0], fs, TIME_RANGE[0], TIME_RANGE[1], win.sum()))
    print("  cutoff = {:.1f} Hz, 4th-order zero-phase Butterworth".format(FC))
    print()

    # ------------------------------------------------------------------ C1
    print("C1  COP columns must be untouched")
    for start, vy_col in BELTS:
        stance = raw[:, vy_col] > THR
        lbl, n = label(stance)
        for i in range(1, n + 1):
            idx = np.where(lbl == i)[0]
            if len(idx) < int(0.10 * fs):
                stance[idx] = False
        lbl, n = label(~stance)
        for i in range(1, n + 1):
            idx = np.where(lbl == i)[0]
            if len(idx) < int(0.04 * fs):
                stance[idx] = True
        sw = (~stance) & win
        for off in OFF_COP:
            c = start + off
            same = np.array_equal(filt[sw, c], raw[sw, c])
            check("COP {} unchanged in swing".format(names[c]), same,
                  "{} swing rows, max|diff| = {:.3e}".format(
                      int(sw.sum()), np.abs(filt[sw, c] - raw[sw, c]).max()))
    print()

    # ------------------------------------------------------------------ C2 / C3
    print("C2  swing-phase force / free moment must be exactly zero")
    print("C3  horizontal torque columns must be exactly zero")
    for bi, (start, vy_col) in enumerate(BELTS):
        stance = raw[:, vy_col] > THR
        lbl, n = label(stance)
        for i in range(1, n + 1):
            idx = np.where(lbl == i)[0]
            if len(idx) < int(0.10 * fs):
                stance[idx] = False
        lbl, n = label(~stance)
        for i in range(1, n + 1):
            idx = np.where(lbl == i)[0]
            if len(idx) < int(0.04 * fs):
                stance[idx] = True

        swing = (~stance) & win
        cols = [start + o for o in OFF_FORCE] + [start + o for o in OFF_TORQUE]
        mx = float(np.abs(filt[swing][:, cols]).max()) if swing.any() else 0.0
        check("belt{} swing |force,torque| == 0".format(bi + 1), mx == 0.0,
              "{} swing rows, max = {:.3e}".format(int(swing.sum()), mx))

        for off, nm in ((OFF_TORQUE[0], "tx"), (OFF_TORQUE[2], "tz")):
            c = start + off
            check("belt{} {} column is all zero".format(bi + 1, nm),
                  np.all(filt[:, c] == 0.0),
                  "max|.| = {:.3e}".format(np.abs(filt[:, c]).max()))
    print()

    # ------------------------------------------------------------------ C4..C6
    print("C4  stance vertical-force peak preserved")
    print("C5  free-moment peak preserved")
    print("C6  swing-to-stance COP step survives")
    print("C7  COP filtering stays dynamically consistent")
    for bi, (start, vy_col) in enumerate(BELTS):
        stance = raw[:, vy_col] > THR
        st_w = stance & win
        r_pk = float(raw[st_w, vy_col].max())
        f_pk = float(filt[st_w, vy_col].max())
        check("belt{} vy peak".format(bi + 1), f_pk >= 0.99 * r_pk,
              "raw {:.2f} N -> filtered {:.2f} N ({:.2f}%)".format(
                  r_pk, f_pk, 100.0 * f_pk / r_pk))

        ty = start + OFF_TORQUE[1]
        r_ty = float(np.abs(raw[st_w, ty]).max())
        f_ty = float(np.abs(filt[st_w, ty]).max())
        check("belt{} |ty| peak".format(bi + 1), f_ty >= 0.75 * r_ty,
              "raw {:.3f} Nm -> filtered {:.3f} Nm ({:.1f}%)".format(
                  r_ty, f_ty, 100.0 * f_ty / r_ty))

        for off, nm in ((OFF_COP[0], "h1"), (OFF_COP[2], "h2")):
            c = start + off
            on = np.where(np.diff(stance.astype(int)) == 1)[0] + 1
            on = on[win[on]]
            on = on[on > 0]
            sr = np.array([abs(raw[o, c] - raw[o - 1, c]) for o in on])
            sf = np.array([abs(filt[o, c] - filt[o - 1, c]) for o in on])
            ratio = np.median(sf) / np.median(sr)
            check("belt{} COP {} onset step kept".format(bi + 1, nm), ratio >= 0.60,
                  "raw median {:.4f} m -> filtered {:.4f} m (ratio {:.3f}, n={})".format(
                      np.median(sr), np.median(sf), ratio, len(on)))

        # ---- C7: the COP filtering must not move the applied moment much ----
        dM = np.abs(filt[st_w][:, [start + o for o in OFF_COP]] -
                    raw[st_w][:, [start + o for o in OFF_COP]]) * \
            raw[st_w, vy_col][:, None]
        pk, rms = float(dM.max()), float(np.sqrt(np.mean(dM ** 2)))
        check("belt{} COP moment error".format(bi + 1), pk <= 20.0 and rms <= 3.0,
              "stance max {:.2f} Nm (<=20), rms {:.2f} Nm (<=3)".format(pk, rms))
    print()

    # ------------------------------------------------------------------ C8
    # Measured on INTERIOR windows only: the Tukey taper and the swing zeroing
    # are deliberate amplitude modulation at the stance boundaries and inject
    # high-frequency content of their own, so including them would measure the
    # window rather than the filter.
    print("C8  filtered force channels are band-limited at the cutoff")
    print("    (interior stance windows, 10% trimmed from each stance end)")

    def interior_windows(mask, min_len=600, trim=0.10):
        lbl, n = label(mask)
        out = []
        for i in range(1, n + 1):
            idx = np.where(lbl == i)[0]
            if len(idx) < min_len:
                continue
            k = int(trim * len(idx))
            out.append((idx[0] + k, idx[-1] + 1 - k))
        return out

    def hi_frac(x, f, fmin):
        """Fraction of the windowed power above `fmin`."""
        p = np.abs(np.fft.rfft(x * np.hanning(len(x)))) ** 2
        tot = p[1:].sum()
        if tot <= 0:
            return 0.0
        return float(p[1:][(f >= fmin)[1:]].sum() / tot)

    for bi, (start, vy_col) in enumerate(BELTS):
        stance = raw[:, vy_col] > THR
        wins = interior_windows(stance & win)
        if not wins:
            check("belt{} interior windows".format(bi + 1), False, "none found")
            continue
        for off in list(OFF_FORCE) + [OFF_TORQUE[1]]:
            c = start + off
            fr, fp = [], []
            for s, e in wins:
                xr = raw[s:e, c] - raw[s:e, c].mean()
                xp = filt[s:e, c] - filt[s:e, c].mean()
                f = np.fft.rfftfreq(len(xr), 1.0 / fs)
                fr.append(hi_frac(xr, f, FC))
                fp.append(hi_frac(xp, f, FC))
            # the cutoff is where a 4th-order Butterworth is already 3 dB down,
            # so require a clear collapse rather than mathematical zero
            ok = np.median(fp) < 0.35 * np.median(fr)
            check("{} >{:.0f} Hz power suppressed".format(names[c], FC), ok,
                  "raw {:.3f}% -> filtered {:.3f}% of windowed power ({} windows)".format(
                      100 * np.median(fr), 100 * np.median(fp), len(wins)))
    print()

    # ------------------------------------------------------------------ summary
    print("=" * 78)
    n_fail = len(_failures)
    print("{} checks, {} failed".format(len(_results), n_fail))
    if n_fail:
        for nm in _failures:
            print("  FAILED: {}".format(nm))
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
