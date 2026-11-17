"""Run a NO-sEMG-tracking MocoInverse baseline over the same 30 s / 4-chunk
structure as the B30_cont2 (sEMG-tracked) run, then compare against it.

Purpose: answer "what changes when we add sEMG tracking?" on identical
settings (model, IK, GRF, mesh, segment/overlap, tol, max_iter, initial guess).

Outputs (all under experiments/moco/outputs/noemg_baseline/):
  logs/nb_u001_c{0..3}.log
  moco_inverse_solution_nb_u001_c{0..3}.sto   (per chunk, merged retained segs)
  moco_inverse_solution_nb_u001_merged.sto    (30 s, stitched)
  noemg_baseline_compare.txt                  (side-by-side metrics)
"""

import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import opensim as osim

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "outputs")
BASE = os.path.join(OUT, "noemg_baseline")
LOGS = os.path.join(BASE, "logs")
os.makedirs(LOGS, exist_ok=True)

PY = sys.executable
CHUNKS = [(210.0, 217.5), (217.5, 225.0), (225.0, 232.5), (232.5, 240.0)]
MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r"]
REF = os.path.join(OUT, "emg_excitation_210_240_right.sto")
TRACKED = os.path.join(OUT, "moco_inverse_solution_B30_cont2_u021_merged.sto")
GUESS = os.path.join(OUT, "moco_inverse_solution_B30_right_final.sto")
TAG = "nb_u001"

MESH = "0.02"
TOL = "1e-3"
MAX_ITER = "3000"


def read_ref(path):
    lines = open(path, encoding="utf-8").readlines()
    ih = [i for i, ln in enumerate(lines)
          if ln.strip().lower().startswith("endheader")][0]
    labels = lines[ih + 1].strip().split("\t")
    data = np.genfromtxt(path, skip_header=ih + 2)
    t = data[:, 0]
    return t, {m: data[:, labels.index("/forceset/" + m)] for m in MUSCLES}


def trapz_weights(t):
    h = t[1] - t[0]
    w = np.ones(len(t)) * h
    w[0] = w[-1] = h / 2.0
    return w


def control(T, m):
    return np.array(T.getControlMat("/forceset/" + m)).reshape(-1)


def metrics(T, ref_t, ref_cols, s):
    """mean cos / pearson / nRMSE of (u vs e) and of (s*u vs e), 5 muscles."""
    t = np.array(T.getTimeMat()).reshape(-1)
    rows_a, rows_b = [], []
    for m, sv in zip(MUSCLES, s):
        u = control(T, m)
        e = np.interp(t, ref_t, ref_cols[m])
        rows_a.append(_corr(u, e))
        rows_b.append(_corr(sv * u, e))
    return np.mean(rows_a, axis=0), np.mean(rows_b, axis=0)


def _corr(u, e):
    c = float(np.dot(u, e) / (np.linalg.norm(u) * np.linalg.norm(e) + 1e-12))
    r = float(np.corrcoef(u, e)[0, 1])
    nr = float(np.sqrt(np.mean((u - e) ** 2)) / np.sqrt(np.mean(e ** 2)))
    return (c, r, nr)


def tracking_obj(T, ref_t, ref_cols, s, weight=20.0):
    t = np.array(T.getTimeMat()).reshape(-1)
    w = trapz_weights(t)
    tot = 0.0
    for m, sv in zip(MUSCLES, s):
        u = control(T, m)
        e = np.interp(t, ref_t, ref_cols[m])
        tot += weight * float(np.sum(w * (e - sv * u) ** 2))
    return tot


def merge(files, out_path):
    trajs = [osim.MocoTrajectory(p) for p in files]
    times, states, controls = [], {}, {}
    for j, T in enumerate(trajs):
        ts = np.array(T.getTimeMat()).reshape(-1)
        times.append(ts if j == 0 else ts[1:])
        for n in T.getStateNames():
            v = np.array([T.getState(n)[i] for i in range(T.getState(n).size())])
            states.setdefault(n, []).append(v if j == 0 else v[1:])
        for n in T.getControlNames():
            v = np.array([T.getControl(n)[i] for i in range(T.getControl(n).size())])
            controls.setdefault(n, []).append(v if j == 0 else v[1:])
    t = np.concatenate(times)
    T = trajs[0]
    T.setNumTimes(len(t))
    T.setTime(t.tolist())
    for n in T.getStateNames():
        T.setState(n, np.concatenate(states[n]).tolist())
    for n in T.getControlNames():
        T.setControl(n, np.concatenate(controls[n]).tolist())
    T.write(out_path)
    return T


def best_s(T, ref_t, ref_cols):
    t = np.array(T.getTimeMat()).reshape(-1)
    w = trapz_weights(t)
    out = []
    for m in MUSCLES:
        u = control(T, m)
        e = np.interp(t, ref_t, ref_cols[m])
        out.append(float(np.sum(w * e * u) / np.sum(w * u * u)))
    return out


def run_chunk(ci):
    start, end = CHUNKS[ci]
    tag = f"{TAG}_c{ci}"
    segdir = os.path.join(BASE, f"segments_{tag}")
    log = os.path.join(LOGS, tag + ".log")
    cmd = [
        PY, os.path.join(HERE, "moco_inverse.py"),
        "--start", str(start), "--end", str(end),
        "--segment", "2.0", "--overlap", "0.5", "--mesh", MESH,
        "--sagittal-cost", "--nonsagittal-weight", "1e-3",
        "--tol", TOL, "--max-iter", MAX_ITER,
        "--emg-initial-guess", GUESS,
        "--seg-dir", segdir,
        "--tag", tag,
    ]
    print("RUN:", " ".join(cmd), flush=True)
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT)
    return ci, p.returncode, time.time() - t0, log, segdir, tag


def main():
    ref_t, ref_cols = read_ref(REF)
    # s from the tracked run (its final outer update)
    ck = json.load(open(os.path.join(OUT, "B30_cont2_checkpoint.json"),
                        encoding="utf-8"))
    s_tracked = [float(x) for x in ck["s_hat"]]

    print("=" * 78, flush=True)
    print("NO-sEMG-tracking baseline for B30 30 s right leg", flush=True)
    print("s (from tracked run) =", s_tracked, flush=True)
    print("guess =", GUESS, flush=True)
    print("=" * 78, flush=True)

    t_all = time.time()
    analyze_only = "--analyze-only" in sys.argv
    merged = os.path.join(BASE, f"moco_inverse_solution_{TAG}_merged.sto")

    if analyze_only:
        if not os.path.isfile(merged):
            print("ERROR: --analyze-only but merged solution missing", flush=True)
            return 1
        T_no = osim.MocoTrajectory(merged)
    else:
        with ThreadPoolExecutor(max_workers=4) as ex:
            results = list(ex.map(run_chunk, range(4)))
        ok = all(rc == 0 for _, rc, _, _, _, _ in results)

        chunk_files = []
        for ci, rc, dt, log, segdir, tag in results:
            print(f"chunk {ci}: rc={rc} wall={dt:.0f}s", flush=True)
            p = os.path.join(OUT, f"moco_inverse_solution_{tag}.sto")
            if os.path.isfile(p):
                chunk_files.append(p)
            else:
                print(f"  WARNING: missing chunk solution {p}", flush=True)

        if not ok or len(chunk_files) != 4:
            print("ERROR: at least one chunk failed or is missing; see logs", flush=True)
            return 1

        # stitch the four retained 7.5 s chunk trajectories into one 30 s trajectory
        T_no = merge(chunk_files, merged)
        print(f"merged -> {merged}", flush=True)
        print(f"merged time range = {T_no.getTimeMat()[0]} .. "
              f"{T_no.getTimeMat()[-1]} ({T_no.getNumTimes()} samples)", flush=True)

    # tracked run, same evaluation function
    T_tr = osim.MocoTrajectory(TRACKED)

    s_best = best_s(T_no, ref_t, ref_cols)
    m_no_raw, m_no_scaled = metrics(T_no, ref_t, ref_cols, s_tracked)
    m_tr_raw, m_tr_scaled = metrics(T_tr, ref_t, ref_cols, s_tracked)
    tr_no = tracking_obj(T_no, ref_t, ref_cols, s_tracked)
    tr_tr = tracking_obj(T_tr, ref_t, ref_cols, s_tracked)

    # effort: sum of excitation_effort over segments, from logs
    def effort_from_logs(tag_prefix):
        eff = 0.0
        for ci in range(4):
            p = os.path.join(LOGS, f"{tag_prefix}_c{ci}.log")
            if not os.path.isfile(p):
                continue
            for ln in open(p, encoding="utf-8", errors="ignore"):
                if ln.strip().startswith("excitation_effort"):
                    eff += float(ln.split(":")[1])
        return eff

    eff_no = effort_from_logs(TAG)

    def effort_from_summaries(prefix):
        tot = 0.0
        for ci in range(4):
            p = os.path.join(OUT, f"moco_inverse_summary_{prefix}_c{ci}.txt")
            if not os.path.isfile(p):
                continue
            for ln in open(p, encoding="utf-8"):
                if ln.startswith("objective="):
                    tot += float(ln.split("=", 1)[1])
                    break
        return tot

    eff_tr = effort_from_summaries("B30_cont2_u021")

    # tracked run: effort / tracking components from its progress log (last row)
    eff_tr_comp = trk_tr_comp = float("nan")
    prog = os.path.join(OUT, "B30_cont2_progress.csv")
    if os.path.isfile(prog):
        rows = list(csv.DictReader(open(prog, encoding="utf-8")))
        if rows:
            eff_tr_comp = float(rows[-1]["effort"])
            trk_tr_comp = float(rows[-1]["tracking"])

    lines = []
    lines.append("# NO-tracking vs tracked MocoInverse baseline (B30 30 s right leg)")
    lines.append(f"generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"tracked solution : {TRACKED}")
    lines.append(f"no-track solution: {merged}")
    lines.append(f"initial guess    : {GUESS}")
    lines.append(f"s (tracked run)  : {s_tracked}")
    lines.append(f"s (best fit, no-track): {['%.4f' % v for v in s_best]}")
    lines.append("")
    lines.append("## same-domain comparison (both stitched 30 s, 3001 samples, 210-240 s)")
    lines.append("metric                        | no-tracking | tracked")
    lines.append("mean cos  (u vs e)            | %11.4f | %7.4f" % (m_no_raw[0], m_tr_raw[0]))
    lines.append("mean Pearson r (u vs e)       | %11.4f | %7.4f" % (m_no_raw[1], m_tr_raw[1]))
    lines.append("mean nRMSE (u vs e)           | %11.4f | %7.4f" % (m_no_raw[2], m_tr_raw[2]))
    lines.append("EMG tracking objective        | %11.2f | %7.2f" % (tr_no, trk_tr_comp))
    lines.append("excitation effort (component) | %11.2f | %7.2f" % (eff_no, eff_tr_comp))
    lines.append("total objective (as reported)  | %10.2f | %7.2f" % (eff_no, eff_tr))
    lines.append("")
    lines.append("## per-muscle (u vs e): cos / r / nRMSE     no-track            tracked")
    for i, m in enumerate(MUSCLES):
        a = _corr(control(T_no, m), np.interp(
            np.array(T_no.getTimeMat()).reshape(-1), ref_t, ref_cols[m]))
        b = _corr(control(T_tr, m), np.interp(
            np.array(T_tr.getTimeMat()).reshape(-1), ref_t, ref_cols[m]))
        lines.append("%-12s no-track %.3f/%.3f/%.3f | tracked %.3f/%.3f/%.3f"
                     % (m, a[0], a[1], a[2], b[0], b[1], b[2]))
    lines.append("")
    lines.append("Note: the tracked run's own log reports mean cos 0.896 / r 0.752 / nRMSE 0.378")
    lines.append("on its 3 s solve-window domain; the numbers above are recomputed on the")
    lines.append("stitched 30 s domain so both columns use identical samples.")
    txt = "\n".join(lines) + "\n"
    out_txt = os.path.join(BASE, "noemg_baseline_compare.txt")
    open(out_txt, "w", encoding="utf-8").write(txt)
    print("\n" + txt, flush=True)
    print(f"[done] total wall = {time.time() - t_all:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
