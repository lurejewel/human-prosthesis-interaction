"""sEMG preprocessing for the JW multi-device EMGworks export.

Python port of experiments/moco/process_semg.m with the same chain and the same
numbers (despike 2.5xP99.9, 10 Hz high-pass, rectify, 5 Hz envelope, per-channel
max normalisation), but the column layout is resolved by header name: the JW
export interleaves force-plate, accelerometer and IMU blocks, so the fixed
12-column textscan of the MATLAB version reads the wrong signals.

Self-contained: with no --data-dir it processes the trial folder it lives in
and writes every output next to itself.  The analysis window defaults to the
one recorded in steady_state_summary.json.
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import butter, filtfilt

CHANNELS = [f"EMG{k}" for k in range(7, 17)]
MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r",
           "rect_fem_l", "vasti_l", "hamstrings_l", "tib_ant_l", "gastroc_l"]

DESPIKE_K = 2.5         # |x| > K * P99.9 flags an artifact
DESPIKE_MARGIN = 10     # samples grown around each flagged run (5 ms)
DESPIKE_NEIGHBORS = 50  # clean support samples per side for the PCHIP rebuild (25 ms)
HP_HZ = 10.0
LP_HZ = 5.0
ORDER = 2
PAD_S = 1.0             # processed beyond the analysis window for solver burn-in


def read_semg(path, wanted):
    """Return (fs, t, columns) reading only the named columns of the export."""
    with open(path, encoding="utf-8-sig") as fh:
        if fh.readline().strip() != "Devices":
            sys.exit(f"{path}: not a Delsys EMGworks 'Devices' export")
        fs = float(fh.readline().strip())
        fh.readline()                       # device block names
        header = fh.readline().rstrip("\n").split(",")
        fh.readline()                       # units row
        idx = {n: i for i, n in enumerate(header)}
        missing = [n for n in wanted if n not in idx]
        if missing:
            sys.exit(f"{path}: missing columns {missing}")
        pick = [idx[n] for n in wanted]
        frames, subs, vals = [], [], []
        for line in fh:
            v = line.split(",")
            if len(v) <= max(pick):
                if frames:
                    break         # blank line ends the 2 kHz block
                continue
            if v[0] == "Frame":
                if frames:
                    break         # a second device block (100 Hz) follows
                continue
            try:
                fr = float(v[0])
                sb = float(v[1])
            except ValueError:
                continue          # units / device rows inside the block
            frames.append(fr)
            subs.append(sb)
            vals.append([float(v[i]) if v[i] else np.nan for i in pick])
    frame = np.asarray(frames)
    sub = np.asarray(subs)
    nsub = int(sub.max()) + 1
    t = ((frame - 1) * nsub + sub) / fs
    if np.any(np.diff(t) <= 0):
        sys.exit(f"{path}: time base not strictly increasing; the export repeats blocks")
    return fs, t, np.asarray(vals)


def despike(x, fs):
    """Replace artifact runs by a PCHIP bridge between clean neighbours."""
    thr = DESPIKE_K * np.percentile(np.abs(x), 99.9)
    bad = np.abs(x) > thr
    if not bad.any():
        return x.copy(), np.zeros(0), thr
    grown = np.convolve(bad.astype(int),
                        np.ones(2 * DESPIKE_MARGIN + 1, int), "same") > 0
    edges = np.flatnonzero(np.diff(np.concatenate(([0], grown.astype(int), [0]))))
    starts, ends = edges[0::2], edges[1::2]
    y = x.copy()
    intervals = []
    for s, e in zip(starts, ends):
        left = np.flatnonzero(~grown[:s])[-DESPIKE_NEIGHBORS:]
        right = np.flatnonzero(~grown[e:])[:DESPIKE_NEIGHBORS] + e
        sup = np.concatenate((left, right))
        if len(sup) >= 4:
            y[s:e] = PchipInterpolator(sup, x[sup])(np.arange(s, e))
        elif len(sup) >= 2:
            y[s:e] = np.interp(np.arange(s, e), sup, x[sup])
        else:
            y[s:e] = 0.0
        intervals.append((s / fs, e / fs, int(e - s), float(np.abs(x[s:e]).max())))
    return y, intervals, thr


def envelope(x, fs):
    b, a = butter(ORDER, HP_HZ / (fs / 2), btype="high")
    y = filtfilt(b, a, x)
    y = np.abs(y)
    b, a = butter(ORDER, LP_HZ / (fs / 2), btype="low")
    y = filtfilt(b, a, y)
    return np.maximum(y, 0.0)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    folder = os.path.abspath(args.data_dir)

    summ = os.path.join(folder, "steady_state_summary.json")
    if args.start is None or args.end is None:
        if not os.path.isfile(summ):
            sys.exit("no --start/--end and no steady_state_summary.json; "
                     "run detect_steady_state.py first")
        w = json.load(open(summ))["analysis_window"]
        args.start = w[0] if args.start is None else args.start
        args.end = w[1] if args.end is None else args.end

    hits = sorted(glob.glob(os.path.join(folder, "*_sEMG.csv")))
    if not hits:
        sys.exit(f"no *_sEMG.csv in {folder}")
    src = hits[0]
    trial = os.path.basename(os.path.normpath(folder))
    stem = os.path.basename(src)[:-len("_sEMG.csv")]

    fs, t, raw = read_semg(src, CHANNELS)
    lo = max(t[0], args.start - PAD_S)
    hi = min(t[-1], args.end + PAD_S)
    sel = (t >= lo) & (t <= hi)
    t, raw = t[sel], raw[sel]
    if np.isnan(raw).any():
        sys.exit("NaN samples inside the analysis window; inspect the export")
    print(f"{trial}: {src}")
    print(f"  fs {fs:.0f} Hz, window [{args.start:.2f}, {args.end:.2f}] s, "
          f"processed [{t[0]:.2f}, {t[-1]:.2f}] s ({len(t)} samples)")

    clean, art, envn = [], [], []
    inner = (t >= args.start) & (t <= args.end)
    for j, ch in enumerate(CHANNELS):
        y, intervals, thr = despike(raw[:, j], fs)
        e = envelope(y, fs)
        scale = e[inner].max()
        if scale <= 0:
            sys.exit(f"{ch}: envelope is identically zero in the window")
        clean.append(y)
        envn.append(e / scale)
        for t0, t1, n, peak in intervals:
            art.append((ch, MUSCLES[j], t0, t1, n, peak, peak / thr))
        print(f"  {ch:5s} -> {MUSCLES[j]:13s} raw|x|max={np.abs(raw[:, j]).max():.5f} V "
              f"despike_thr={thr:.5f} V runs={len(intervals)} envmax={scale:.5f} V")
    envn = np.column_stack(envn)

    out_csv = os.path.join(folder, f"{stem}_semg_processed.csv")
    with open(out_csv, "w") as fh:
        fh.write("Time_s," + ",".join(CHANNELS) + "\n")
        for i in range(len(t)):
            fh.write(f"{t[i]:.6f}," + ",".join(f"{v:.8f}" for v in envn[i]) + "\n")
    print(f"  wrote {out_csv}")

    rep = os.path.join(folder, f"{stem}_semg_artifacts.csv")
    with open(rep, "w") as fh:
        fh.write("EMG,Muscle,t_start_s,t_end_s,SamplesFlagged,PeakRemoved_V,PeakOverThr\n")
        for r in art:
            fh.write(f"{r[0]},{r[1]},{r[2]:.4f},{r[3]:.4f},{r[4]},{r[5]:.6f},{r[6]:.2f}\n")
    print(f"  wrote {rep}  ({len(art)} artifact runs)")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(5, 2, figsize=(14, 10), sharex=True)
        for j, ax in enumerate(axes.T.ravel()):
            ax.plot(t, envn[:, j], lw=0.6)
            ax.set_title(f"{CHANNELS[j]} {MUSCLES[j]}", fontsize=9)
            ax.set_ylim(0, 1.05)
            ax.axvspan(args.start, args.end, color="tab:blue", alpha=0.08)
        for ax in axes[-1]:
            ax.set_xlabel("time (s)")
        fig.suptitle(f"{trial} normalised sEMG envelope (10 Hz HP / rectify / 5 Hz LP)", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(folder, f"{stem}_semg_envelope.png"), dpi=130)
        plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())