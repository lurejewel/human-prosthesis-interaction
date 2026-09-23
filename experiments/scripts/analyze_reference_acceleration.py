#!/usr/bin/env python3
"""Quantify the high-frequency content of the reference acceleration fed to CMC.

CMC never differentiates the measured kinematics.  CMCTool.cpp builds a
GCVSplineSet (degree 5) through the desired kinematics and takes the desired
acceleration from uDotSet = d/dt(uSet), i.e. the *second derivative of an
interpolating spline*.  With `use_fast_optimization_target=true` the resulting
qdd* enters the static optimisation as an exact equality constraint, so any
high-frequency content in qdd* is propagated straight into the muscle forces.

This script answers, on real data:

  1. How much high-frequency energy does the spline acceleration carry,
     compared with a properly smoothed reference?
  2. Does the acceleration stay consistent with the position it came from?
  3. Where is the muscle activation bandwidth, so where should the cut be?

Usage
-----
    python analyze_reference_acceleration.py                 # both files
    python analyze_reference_acceleration.py --plot          # + figures

Run it with an interpreter that has OpenSim, numpy and scipy, e.g.
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# numpy 2.0 renamed trapz -> trapezoid; keep working on both
_trapz = getattr(np, "trapezoid", None) or np.trapz

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

DEFAULT_FILES = [
    os.path.join(DATA_DIR, "level_walking_ik_filtered.mot"),  # 15 Hz Butterworth
    os.path.join(DATA_DIR, "level_walking_ik.mot"),           # raw
]

# Muscle activation dynamics: the fastest mode a muscle can follow.
TAU_ACT = 0.015     # s, activation time constant  -> ~10.6 Hz corner
TAU_DEACT = 0.050   # s, deactivation time constant -> ~3.2 Hz corner

# Reference trajectory band of interest: walking has almost no content above
# 6 Hz (Winter; ~5 harmonics of ~1 Hz), so anything above this in qdd* is
# numerical, not physiological.
SIGNAL_BAND_HZ = 6.0

ANALYSIS_BANDS_HZ = [10.0, 20.0, 30.0]


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_mot(path):
    """Read an OpenSim .mot/.sto file.  Returns (labels, time, data, in_degrees).

    `data` is (n_samples, n_labels); the time column is returned separately.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    in_degrees = False
    i = 0
    for i, line in enumerate(lines):
        low = line.strip().lower()
        if low.startswith("indegrees"):
            in_degrees = "yes" in low
        if low.startswith("endheader"):
            break
    else:
        sys.exit("no 'endheader' found in {}".format(path))

    labels = lines[i + 1].split()
    rows = [ln.split() for ln in lines[i + 2:] if ln.strip()]
    arr = np.array(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != len(labels):
        sys.exit("ragged data in {}: {} labels, {} columns"
                 .format(path, len(labels), arr.shape[1] if arr.ndim == 2 else "?"))
    return labels[1:], arr[:, 0], arr[:, 1:], in_degrees


# ---------------------------------------------------------------------------
# Reference construction -- mirrors CMCTool's GCVSplineSet(5, store)
# ---------------------------------------------------------------------------
def spline_reference(t, x, degree=5):
    """Interpolating spline through every sample, as OpenSim's GCVSplineSet
    does for dense data, returned with its analytic derivatives."""
    from scipy.interpolate import make_interp_spline
    k = min(degree, len(t) - 1)
    spl = make_interp_spline(t, x, k=k)
    dspl = spl.derivative(1)
    ddspl = spl.derivative(2)
    return spl(t), dspl(t), ddspl(t)


def finite_difference_accel(t, x):
    """Second-order accurate, non-uniform finite-difference acceleration."""
    return np.gradient(np.gradient(x, t, axis=0), t, axis=0)


def smooth_reference(t, x, cutoff_hz):
    """Zero-phase Butterworth reference.

    CMC's own filter (`lowpass_cutoff_frequency`, FIR order 50) is applied to
    the *positions* before the spline is built.  That is the right place: the
    acceleration is then the derivative of an already smooth trajectory, so
    q, qd and qdd stay mutually consistent.
    """
    from scipy.signal import butter, filtfilt
    fs = 1.0 / np.mean(np.diff(t))
    nyq = 0.5 * fs
    if cutoff_hz <= 0 or cutoff_hz >= nyq:
        return x, None
    b, a = butter(4, cutoff_hz / nyq, btype="low")
    padlen = min(3 * max(len(a), len(b)), len(x) - 1)
    # axis=0: filter along time, not across coordinates
    xs = filtfilt(b, a, x, axis=0, padlen=padlen)
    return xs, (b, a)


# ---------------------------------------------------------------------------
# Spectral metrics
# ---------------------------------------------------------------------------
def spectrum(t, x):
    from scipy.signal import welch
    fs = 1.0 / np.mean(np.diff(t))
    nper = min(len(x), 2048)
    f, pxx = welch(x, fs=fs, nperseg=nper)
    return f, pxx


def band_energy_fraction(f, pxx, f_lo, f_hi=None):
    """Fraction of total power above `f_lo` (or inside [f_lo, f_hi])."""
    total = _trapz(pxx, f)
    if total <= 0:
        return 0.0
    if f_hi is None:
        m = f >= f_lo
    else:
        m = (f >= f_lo) & (f <= f_hi)
    return float(_trapz(pxx[m], f[m]) / total)


def consistency_error(t, q, qd, qdd):
    """Reconstruct q by integrating qdd twice; report peak discrepancy.

    A large value means qdd is NOT the second derivative of the q that the
    controller is trying to track -- i.e. the PD demand is fighting itself.
    """
    qd_rec = qd[0] + np.concatenate(([0.0], np.cumsum(0.5 * (qdd[1:] + qdd[:-1])
                                                    * np.diff(t))))
    q_rec = q[0] + np.concatenate(([0.0], np.cumsum(0.5 * (qd_rec[1:] + qd_rec[:-1])
                                                   * np.diff(t))))
    # remove the linear drift that any double integration accumulates
    A = np.vstack([t - t[0], np.ones_like(t)]).T
    coef, *_ = np.linalg.lstsq(A, q_rec - q, rcond=None)
    resid = (q_rec - q) - A @ coef
    return float(np.max(np.abs(resid))), float(np.sqrt(np.mean(resid ** 2)))


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def analyze(path, cutoff=SIGNAL_BAND_HZ, do_plot=False):
    labels, t, data, in_degrees = read_mot(path)
    unit = "deg" if in_degrees else "rad"
    fs = 1.0 / np.mean(np.diff(t))
    print("=" * 78)
    print(os.path.basename(path))
    print("=" * 78)
    print("  samples          : {}  (t = {:.2f} .. {:.2f} s, fs = {:.0f} Hz)"
          .format(len(t), t[0], t[-1], fs))
    print("  coordinates      : {}  [{}]".format(len(labels), unit))
    print("  CMC control rate : {:.0f} Hz  (cmc_time_window = 0.010 s)"
          .format(1.0 / 0.010))
    print("  muscle bandwidth : activation {:.1f} Hz, deactivation {:.1f} Hz"
          .format(1.0 / (2 * np.pi * TAU_ACT), 1.0 / (2 * np.pi * TAU_DEACT)))
    print()

    # ---- what CMC actually feeds into the PD law -------------------------
    _, qd_s, qdd_s = spline_reference(t, data)
    qdd_fd = finite_difference_accel(t, data)

    # ---- the same reference, smoothed before splining --------------------
    data_s, _ = smooth_reference(t, data, cutoff)
    _, qd_sm, qdd_sm = spline_reference(t, data_s)

    header = ("  {:<22} {:>12} {:>12} {:>12}"
              .format("coordinate", "rms qdd", "rms accel", ">20 Hz"))
    print(header)
    print("  " + "-" * (len(header) - 2))
    print("  {:<22} {:>12} {:>12} {:>12}"
          .format("", "spline", "finite diff", "power frac"))
    print("  " + "-" * (len(header) - 2))

    rows = []
    for j, name in enumerate(labels):
        qdd_fd_col = qdd_fd[:, j]
        # compare only where the finite difference is meaningful
        good = np.isfinite(qdd_fd_col)
        rms_spline = float(np.sqrt(np.mean(qdd_s[:, j] ** 2)))
        rms_fd = float(np.sqrt(np.mean(qdd_fd_col[good] ** 2)))
        f, pxx = spectrum(t, qdd_s[:, j])
        frac20 = band_energy_fraction(f, pxx, 20.0)
        rows.append((name, rms_spline, rms_fd, frac20))

    # print the noisiest coordinates first -- those are the ones that decide
    # whether the QP target is achievable
    for name, rms_spline, rms_fd, frac20 in sorted(rows, key=lambda r: -r[3]):
        flag = ""
        if frac20 > 0.05:
            flag = "  <-- noisy target"
        print("  {:<22} {:>12.4g} {:>12.4g} {:>11.1f}%{}"
              .format(name, rms_spline, rms_fd, 100 * frac20, flag))

    print()
    print("  Aggregate high-frequency power fraction of the spline acceleration:")
    for lo in ANALYSIS_BANDS_HZ:
        fr = []
        for j in range(len(labels)):
            f, pxx = spectrum(t, qdd_s[:, j])
            fr.append(band_energy_fraction(f, pxx, lo))
        print("    above {:>4.0f} Hz : mean {:>6.2f}%   max {:>6.2f}%"
              .format(lo, 100 * np.mean(fr), 100 * np.max(fr)))

    # ---- the demand the muscles physically CANNOT follow ---------------
    # A muscle cannot change force on a timescale faster than its activation
    # time constant, while CMC demands a *new constant* excitation every
    # control window.  Compare the two directly: the reference acceleration has
    # an rms value A over a band, its rate of change has rms A_dot, so within
    # one control window the demanded change is A_dot*dt; a first-order lag
    # with time constant tau can deliver at most (1-exp(-dt/tau)) of a step.
    dt_win = 0.010                      # cmc_time_window
    alpha = {}
    for label, tau in (("activation", TAU_ACT), ("deactivation", TAU_DEACT)):
        alpha[label] = 1.0 - np.exp(-dt_win / tau)

    print()
    print("  Can the muscles deliver the demanded acceleration change?")
    print("    control window dt = {:.0f} ms".format(1000 * dt_win))
    print("    activation   tau = {:>4.0f} ms -> can cover {:>4.1f}% of a step"
          .format(1000 * TAU_ACT, 100 * alpha["activation"]))
    print("    deactivation tau = {:>4.0f} ms -> can cover {:>4.1f}% of a step"
          .format(1000 * TAU_DEACT, 100 * alpha["deactivation"]))
    print()
    print("    {:<22} {:>12} {:>14} {:>14}".format(
        "coordinate", "rms qdd", "demanded/window", "deliverable"))
    print("    " + "-" * 66)
    req_all, del_all = [], []
    table = []
    for j, name in enumerate(labels):
        qdd = qdd_s[:, j]
        qddd = np.gradient(qdd, t, axis=0)
        a_rms = float(np.sqrt(np.mean(qdd ** 2)))
        adot_rms = float(np.sqrt(np.mean(qddd ** 2)))
        req = adot_rms * dt_win                       # demanded change
        # the slowest (deactivation) mode bounds what can actually be tracked
        deliverable = a_rms * alpha["deactivation"]
        req_all.append(req)
        del_all.append(deliverable)
        table.append((name, a_rms, req, deliverable))
    for name, a_rms, req, deliverable in sorted(table, key=lambda r: -r[2] / r[1]):
        ratio = req / a_rms if a_rms > 0 else 0.0
        mark = "  <-- beyond bandwidth" if req > deliverable else ""
        print("    {:<22} {:>12.4g} {:>14.4g} {:>14.4g}  ({:.0f}% of rms){}"
              .format(name, a_rms, req, deliverable, 100 * ratio, mark))
    print("    {:<22} {:>12} {:>14.4g} {:>14.4g}".format(
        "MEAN", "", float(np.mean(req_all)), float(np.mean(del_all))))

    # ---- how low must the reference bandwidth go? ----------------------
    # Absolute rms left in the acceleration above the muscle corner is the
    # quantity that matters: that is how much of qdd* the QP keeps chasing
    # without ever achieving it.
    f_corner = 1.0 / (2 * np.pi * TAU_ACT)
    print()
    print("  Reference bandwidth sweep (absolute rms acceleration above "
          "{:.1f} Hz):".format(f_corner))
    print("    {:<8} {:>14} {:>14} {:>16} {:>12}".format(
        "cutoff", "rms >corner", "jerk rms", "demanded/window", "max |dq|"))
    for fc in (20.0, 15.0, 10.0, 6.0, 4.0, 3.0):
        dat_f, _ = smooth_reference(t, data, fc)
        _, _, qdd_f = spline_reference(t, dat_f)
        hi, jk = [], []
        for j in range(len(labels)):
            f, pxx = spectrum(t, qdd_f[:, j])
            m = f >= f_corner
            hi.append(np.sqrt(max(_trapz(pxx[m], f[m]), 0.0)))
            jk.append(float(np.sqrt(np.mean(np.gradient(qdd_f[:, j], t, axis=0) ** 2))))
        print("    {:<8.0f} {:>13.4g} {:>14.4g} {:>16.4g} {:>12.4g}".format(
            fc, float(np.mean(hi)), float(np.mean(jk)),
            float(np.mean(jk)) * dt_win, np.max(np.abs(dat_f - data))))

    print()
    print("  Consistency of the acceleration reference (double integration residue):")
    worst = 0.0
    worst_name = ""
    for j, name in enumerate(labels):
        pk, rms = consistency_error(t, data[:, j], qd_s[:, j], qdd_s[:, j])
        if pk > worst:
            worst, worst_name = pk, name
    print("    worst peak residue : {:.4g} {}  ({})".format(worst, unit, worst_name))

    # how much does smoothing change the reference?  (tracking vs smoothness)
    err = data_s - data
    print("    effect of {:.0f} Hz smoothing : max |dq| = {:.4g} {}, rms = {:.4g} {}"
          .format(cutoff, np.max(np.abs(err)), unit,
                  np.sqrt(np.mean(err ** 2)), unit))

    # ---- how low must the reference bandwidth go? ----------------------
    # If a muscle cannot follow a mode, asking for it only buys excitation
    # chatter.  Sweep the reference cutoff and report how much of qdd* is
    # left above the muscle corner.
    print()
    print("  Reference bandwidth sweep (fraction of qdd power above 10.6 Hz):")
    print("    {:<10} {:>12} {:>12} {:>14}"
          .format("cutoff", "mean >10.6Hz", "max >10.6Hz", "max |dq| ({})".format(unit)))
    for fc in (20.0, 15.0, 10.0, 6.0, 4.0, 3.0):
        dat_f, _ = smooth_reference(t, data, fc)
        _, _, qdd_f = spline_reference(t, dat_f)
        fr = []
        for j in range(len(labels)):
            f, pxx = spectrum(t, qdd_f[:, j])
            fr.append(band_energy_fraction(f, pxx, 1.0 / (2 * np.pi * TAU_ACT)))
        print("    {:<10.0f} {:>11.2f}% {:>11.2f}% {:>14.4g}"
              .format(fc, 100 * np.mean(fr), 100 * np.max(fr),
                      np.max(np.abs(dat_f - data))))

    print()
    print("  Same statistics for the smoothed reference (what the PD law would see):")
    for lo in ANALYSIS_BANDS_HZ:
        fr = []
        for j in range(len(labels)):
            f, pxx = spectrum(t, qdd_sm[:, j])
            fr.append(band_energy_fraction(f, pxx, lo))
        print("    above {:>4.0f} Hz : mean {:>6.2f}%   max {:>6.2f}%"
              .format(lo, 100 * np.mean(fr), 100 * np.max(fr)))

    if do_plot:
        plot(path, t, labels, qdd_s, qdd_sm, cutoff)

    return {"path": path, "fs": fs, "labels": labels, "qdd_spline": qdd_s,
            "qdd_smooth": qdd_sm}


def plot(path, t, labels, qdd_s, qdd_sm, cutoff):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print("  (matplotlib unavailable: {})".format(exc))
        return

    fs = 1.0 / np.mean(np.diff(t))
    # pick the three coordinates with the most high-frequency content
    scores = []
    for j in range(len(labels)):
        f, pxx = spectrum(t, qdd_s[:, j])
        scores.append((band_energy_fraction(f, pxx, 20.0), j))
    picks = [j for _, j in sorted(scores, reverse=True)[:3]]

    fig, axes = plt.subplots(len(picks), 2, figsize=(13, 3.2 * len(picks)))
    axes = np.atleast_2d(axes)
    for r, j in enumerate(picks):
        ax = axes[r, 0]
        ax.plot(t, qdd_s[:, j], lw=0.6, label="spline qdd (fed to CMC)")
        ax.plot(t, qdd_sm[:, j], lw=1.0,
                label="smoothed {:.0f} Hz then splined".format(cutoff))
        ax.set_title("{}  acceleration".format(labels[j]), fontsize=9)
        ax.set_xlabel("time (s)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)

        ax = axes[r, 1]
        f1, p1 = spectrum(t, qdd_s[:, j])
        f2, p2 = spectrum(t, qdd_sm[:, j])
        ax.loglog(f1[1:], p1[1:], lw=0.8, label="spline")
        ax.loglog(f2[1:], p2[1:], lw=0.8, label="smoothed")
        for fc in (1.0 / (2 * np.pi * TAU_DEACT), 1.0 / (2 * np.pi * TAU_ACT),
                   0.5 * fs):
            ax.axvline(fc, color="k", ls=":", lw=0.7)
        ax.set_title("{}  spectrum".format(labels[j]), fontsize=9)
        ax.set_xlabel("frequency (Hz)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    stem = os.path.splitext(os.path.basename(path))[0]
    out = os.path.join(os.path.dirname(path), stem + "_reference_accel.png")
    fig.savefig(out, dpi=130)
    print("  wrote {}".format(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", default=None,
                    help="kinematics .mot/.sto to analyse (default: both IK files)")
    ap.add_argument("--cutoff", type=float, default=SIGNAL_BAND_HZ,
                    help="cutoff (Hz) for the smoothed comparison reference")
    ap.add_argument("--plot", action="store_true", help="write spectrum figures")
    args = ap.parse_args()

    files = args.files or [f for f in DEFAULT_FILES if os.path.exists(f)]
    if not files:
        sys.exit("no input files found under {}".format(DATA_DIR))
    for p in files:
        if not os.path.exists(p):
            sys.exit("input file not found: {}".format(p))
        analyze(p, cutoff=args.cutoff, do_plot=args.plot)


if __name__ == "__main__":
    main()
