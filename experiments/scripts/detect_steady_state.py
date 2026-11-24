"""Locate the steady-state portion of a treadmill gait trial.

Stability is judged from the GRF alone (heel-strike timing regularity) and
cross-checked against pelvis drift from inverse kinematics.  The analysis
window handed to MocoInverse is the last ``--window`` seconds of the longest
stable span.

The script is self-contained: with no ``--data-dir`` it operates on the
directory it lives in, so a copy placed inside a trial folder reads that
trial's ``*_grf.mot`` / ``*_ik.mot`` and writes every output next to itself.
"""

import argparse
import glob
import json
import os
import sys
from datetime import date

import numpy as np

FORCE_THR = 20.0        # N, vertical force marking ground contact
MIN_STANCE = 0.10       # s, shorter contact runs are noise
MIN_GAP = 0.04          # s, shorter swing gaps are noise
WIN = 10.0              # s, sliding evaluation window
HOP = 1.0               # s, sliding step
CV_MAX = 0.05           # max coefficient of variation of step intervals
RATE_TOL = 0.10         # max relative deviation of local step rate from trial median
DRIFT_MAX = 0.25        # m/s, max pelvis translation drift on the belt
MIN_STEPS = 8           # min steps inside an evaluation window


def read_mot(path):
    with open(path) as fh:
        lines = fh.readlines()
    end = next(i for i, l in enumerate(lines) if l.strip() == "endheader")
    names = lines[end + 1].split()
    data = np.loadtxt(lines[end + 2:])
    return names, data


def read_ik(path):
    names, data = read_mot(path)
    return names, data


def contact_runs(mask, fs):
    """Debounce a contact mask and return (start, end) index pairs."""
    idx = np.flatnonzero(np.diff(np.concatenate(([0], mask.astype(int), [0]))))
    starts, ends = idx[0::2], idx[1::2]
    min_len = int(round(MIN_STANCE * fs))
    keep = (ends - starts) >= min_len
    starts, ends = starts[keep], ends[keep]
    # bridge short swing gaps
    merged_s, merged_e = [], []
    for s, e in zip(starts, ends):
        if merged_s and (s - merged_e[-1]) < int(round(MIN_GAP * fs)):
            merged_e[-1] = e
        else:
            merged_s.append(s)
            merged_e.append(e)
    return list(zip(merged_s, merged_e))


def belt_events(t, vy, fs):
    runs = contact_runs(vy > FORCE_THR, fs)
    # a run touching either record edge is quiet standing, not a step
    runs = [(s, e) for s, e in runs if s > 0 and e < len(t)]
    onsets = np.array([t[s] for s, e in runs])
    durations = np.array([(e - s) / fs for s, e in runs])
    return onsets, durations, runs


def sliding_metrics(t_step, pelvis_t, pelvis_x, t_end):
    """Per-window step-timing regularity and pelvis drift."""
    if len(t_step) < 3:
        return []
    intervals = np.diff(t_step)
    t_mid = 0.5 * (t_step[1:] + t_step[:-1])
    global_median = float(np.median(intervals))
    rows = []
    w0 = t_step[0]
    while w0 + WIN <= t_end + 1e-9:
        w1 = w0 + WIN
        sel = (t_mid >= w0) & (t_mid < w1)
        iv = intervals[sel]
        if len(iv) >= MIN_STEPS - 1:
            mean_iv = float(iv.mean())
            cv = float(iv.std(ddof=1) / mean_iv)
            rate_dev = abs(mean_iv - global_median) / global_median
            pm = (pelvis_t >= w0) & (pelvis_t < w1)
            drift = float(np.polyfit(pelvis_t[pm], pelvis_x[pm], 1)[0]) if pm.sum() > 5 else np.nan
            stable = (cv <= CV_MAX) and (rate_dev <= RATE_TOL) and (abs(drift) <= DRIFT_MAX)
            rows.append(dict(t0=w0, t1=w1, n_steps=int(len(iv) + 1), mean_step=mean_iv,
                             cv=cv, rate_dev=rate_dev, drift=drift, stable=bool(stable)))
        w0 += HOP
    return rows


def longest_stable_span(rows):
    best, cur = None, None
    for r in rows:
        if r["stable"]:
            cur = [r["t0"], r["t1"]] if cur is None else [cur[0], r["t1"]]
            if best is None or (cur[1] - cur[0]) > (best[1] - best[0]):
                best = list(cur)
        else:
            cur = None
    return best


def find_one(pattern, folder):
    hits = sorted(glob.glob(os.path.join(folder, pattern)))
    if not hits:
        sys.exit(f"no file matching {pattern} in {folder}")
    return hits[0]


def analyse(folder, window_len, write_notes):
    grf_path = find_one("*_grf.mot", folder)
    ik_path = find_one("*_ik.mot", folder)
    gnames, g = read_mot(grf_path)
    inames, ik = read_ik(ik_path)
    t = g[:, 0]
    fs = float(1.0 / np.median(np.diff(t)))
    # belt 1 = left, belt 2 = right (see the trial's ExternalLoads XML)
    l_on, l_dur, _ = belt_events(t, g[:, gnames.index("ground_force1_vy")], fs)
    r_on, r_dur, _ = belt_events(t, g[:, gnames.index("ground_force2_vy")], fs)
    steps = np.sort(np.concatenate((l_on, r_on)))

    px = ik[:, inames.index("pelvis_tx")]
    rows = sliding_metrics(steps, ik[:, 0], px, t[-1])
    span = longest_stable_span(rows)
    if span is None:
        sys.exit(f"{folder}: no stable span found; loosen CV_MAX / RATE_TOL")
    stable_len = span[1] - span[0]
    if stable_len < window_len:
        print(f"WARNING: stable span {stable_len:.1f} s < requested {window_len:.1f} s", file=sys.stderr)
    a1 = span[1]
    a0 = max(span[0], a1 - window_len)

    trial_median = float(np.median(np.diff(steps)))
    in_win = (steps >= a0) & (steps <= a1)
    summary = dict(
        trial=os.path.basename(os.path.normpath(folder)),
        grf_file=os.path.basename(grf_path), ik_file=os.path.basename(ik_path),
        grf_fs=fs, ik_fs=float(1.0 / np.median(np.diff(ik[:, 0]))),
        trial_span=[float(t[0]), float(t[-1])],
        n_steps_left=int(len(l_on)), n_steps_right=int(len(r_on)),
        median_step_interval=trial_median,
        cadence_steps_per_min=60.0 / trial_median,
        mean_stance_left=float(l_dur.mean()) if len(l_dur) else None,
        mean_stance_right=float(r_dur.mean()) if len(r_dur) else None,
        stable_span=span, stable_duration=stable_len,
        analysis_window=[a0, a1], analysis_duration=a1 - a0,
        n_windows=len(rows), n_stable_windows=sum(r["stable"] for r in rows),
        params=dict(force_thr=FORCE_THR, win=WIN, hop=HOP, cv_max=CV_MAX,
                    rate_tol=RATE_TOL, drift_max=DRIFT_MAX),
    )

    base = os.path.join(folder, "steady_state")
    with open(base + "_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    if rows:
        keys = list(rows[0].keys())
        with open(base + "_windows.csv", "w") as fh:
            fh.write(",".join(keys) + "\n")
            for r in rows:
                fh.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    np.savetxt(base + "_steps.csv",
               np.column_stack((steps[:-1], np.diff(steps))),
               delimiter=",", header="step_time_s,step_interval_s", comments="")
    plot(folder, t, g, gnames, steps, rows, span, (a0, a1), ik, inames, base + ".png")
    if write_notes:
        write_notes_file(folder, summary)
    return summary


def plot(folder, t, g, gnames, steps, rows, span, win, ik, inames, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    ax[0].plot(t, g[:, gnames.index("ground_force1_vy")], lw=0.4, label="left belt")
    ax[0].plot(t, g[:, gnames.index("ground_force2_vy")], lw=0.4, label="right belt")
    ax[0].set_ylabel("vertical GRF (N)")
    ax[0].legend(loc="upper right", fontsize=8)
    for a in ax:
        a.axvspan(span[0], span[1], color="tab:green", alpha=0.12)
        a.axvspan(win[0], win[1], color="tab:blue", alpha=0.18)
    iv = np.diff(steps)
    ax[1].plot(steps[:-1], iv, ".-", ms=3, lw=0.6)
    ax[1].axhline(np.median(iv), color="k", ls="--", lw=0.8, label="trial median")
    ax[1].set_ylabel("step interval (s)")
    ax[1].legend(fontsize=8)
    if rows:
        wt = [0.5 * (r["t0"] + r["t1"]) for r in rows]
        ax[2].plot(wt, [r["cv"] for r in rows], label="step-interval CV")
        ax[2].plot(wt, [r["rate_dev"] for r in rows], label="rate deviation")
        ax[2].axhline(CV_MAX, color="r", ls="--", lw=0.8)
        ax2 = ax[2].twinx()
        ax2.plot(wt, [r["drift"] for r in rows], color="tab:purple", lw=0.8, label="pelvis drift")
        ax2.set_ylabel("pelvis drift (m/s)", color="tab:purple")
        ax[2].set_ylabel("fraction")
        ax[2].legend(loc="upper left", fontsize=8)
    ax[2].set_xlabel("time (s)")
    fig.suptitle(f"{os.path.basename(os.path.normpath(folder))}  stable [{span[0]:.2f}, {span[1]:.2f}] s"
                 f"   analysis [{win[0]:.2f}, {win[1]:.2f}] s", fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def write_notes_file(folder, s):
    path = os.path.join(folder, "notes.txt")
    old = ""
    if os.path.isfile(path):
        old = open(path).read()
    head = (f"range of stable walking:[{s['stable_span'][0]:.2f}, {s['stable_span'][1]:.2f}] "
            f"(choose the last 30 seconds of stable movement as analysis window)")
    body = f"""
MocoInverse analysis window:[{s['analysis_window'][0]:.2f}, {s['analysis_window'][1]:.2f}]  ({s['analysis_duration']:.2f} s)

[steady-state detection]  generated {date.today().isoformat()} by detect_steady_state.py
method              : heel-strike timing regularity from GRF, {WIN:.0f} s sliding window (hop {HOP:.0f} s)
stable criteria     : step-interval CV <= {CV_MAX:.2f}, local step rate within {RATE_TOL:.0%} of trial median,
                      |pelvis_tx drift| <= {DRIFT_MAX:.2f} m/s, >= {MIN_STEPS} steps per window
contact threshold   : {FORCE_THR:.0f} N vertical, stance >= {MIN_STANCE:.2f} s, swing gap >= {MIN_GAP:.2f} s
trial span          : {s['trial_span'][0]:.2f} - {s['trial_span'][1]:.2f} s  (GRF {s['grf_fs']:.0f} Hz, IK {s['ik_fs']:.0f} Hz)
steps detected      : left {s['n_steps_left']}, right {s['n_steps_right']}
median step interval: {s['median_step_interval']:.4f} s  (cadence {s['cadence_steps_per_min']:.1f} steps/min)
mean stance duration: left {s['mean_stance_left']:.3f} s, right {s['mean_stance_right']:.3f} s
stable span         : {s['stable_span'][0]:.2f} - {s['stable_span'][1]:.2f} s  ({s['stable_duration']:.2f} s, {s['n_stable_windows']}/{s['n_windows']} windows stable)
files               : steady_state_summary.json, steady_state_windows.csv, steady_state_steps.csv, steady_state.png
"""
    with open(path, "w") as fh:
        fh.write(head + "\n" + body)
    return old


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--window", type=float, default=30.0, help="analysis window length (s)")
    ap.add_argument("--no-notes", action="store_true", help="do not rewrite notes.txt")
    ap.add_argument("--all", action="store_true", help="treat --data-dir as a parent of trial folders")
    args = ap.parse_args()

    if args.all:
        folders = sorted(d for d in glob.glob(os.path.join(args.data_dir, "*"))
                         if os.path.isdir(d) and glob.glob(os.path.join(d, "*_grf.mot")))
    else:
        folders = [args.data_dir]
    for f in folders:
        s = analyse(f, args.window, not args.no_notes)
        print(f"{s['trial']:26s} span {s['trial_span'][1]:7.2f}s  "
              f"stable [{s['stable_span'][0]:7.2f},{s['stable_span'][1]:7.2f}] "
              f"({s['stable_duration']:6.2f}s)  analysis [{s['analysis_window'][0]:7.2f},"
              f"{s['analysis_window'][1]:7.2f}]  cadence {s['cadence_steps_per_min']:5.1f}/min")


if __name__ == "__main__":
    main()