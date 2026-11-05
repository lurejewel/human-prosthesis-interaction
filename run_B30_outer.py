"""Continue the 30 s right-leg B outer loop from round 3 to round 13.

Each round:
  - solve 4 chunks in parallel with the current 5 per-muscle s;
  - merge the 4 chunk solutions to a 30 s trajectory;
  - update s_i = sum(e_i u_i) / sum(u_i^2) from all retained samples;
  - append J/tracking/R_s/r_J to B30_right_outer_progress.csv.

Run with the opensim_scripting Python.
"""

import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import opensim as osim

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "experiments", "moco", "outputs")
PY = sys.executable
MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r"]
CHUNKS = [(210.0, 217.5), (217.5, 225.0), (225.0, 232.5), (232.5, 240.0)]
REF = os.path.join(OUT, "emg_excitation_210_240_right.sto")

S = [1.453894, 1.567464, 0.790150, 1.149535, 0.930330]
GUESS = os.path.join(OUT, "moco_inverse_solution_B30_right_final.sto")
PREV_J = 204.111790
START_ROUND, END_ROUND = 3, 13

# reference
lines = open(REF, encoding="utf-8").readlines()
ih = [k for k, line in enumerate(lines)
      if line.strip().lower().startswith("endheader")][0]
labels = lines[ih + 1].strip().split("\t")
ref = np.genfromtxt(REF, skip_header=ih + 2)
ref_t = ref[:, 0]
ref_col = {m: ref[:, labels.index("/forceset/" + m)] for m in MUSCLES}

prog = os.path.join(OUT, "B30_right_outer_progress.csv")
if not os.path.isfile(prog):
    with open(prog, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerow([
            "round", "J", "mean_cos", "mean_pearson", "mean_nRMSE",
            "s_rect", "s_vasti", "s_ham", "s_tib", "s_gastroc",
            "R_s", "rJ",
        ])


def solve_chunk(round_no, ci, s_vals, guess):
    tag = f"B30_right_r{round_no}_chunk{ci}"
    start, end = CHUNKS[ci]
    cmd = [
        PY, os.path.join(ROOT, "experiments", "moco", "moco_inverse.py"),
        "--start", str(start), "--end", str(end),
        "--segment", "2.0", "--overlap", "0.5", "--mesh", "0.02",
        "--sagittal-cost", "--nonsagittal-weight", "1e-3",
        "--emg-tracking", REF, "--emg-weight", "20", "--emg-leg", "right",
        "--emg-native-b-fixed",
        "--emg-s-values", ",".join(f"{v:.6f}" for v in s_vals),
        "--emg-initial-guess", guess,
        "--emg-sparsity-random", "--tol", "1e-4", "--max-iter", "5000",
        "--tag", tag,
    ]
    log = os.path.join(ROOT, f"tmp_{tag}.log")
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run(cmd, cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT)
    return tag, p.returncode


def merge(tags, out_path):
    trajs = [osim.MocoTrajectory(os.path.join(OUT, f"moco_inverse_solution_{t}.sto"))
             for t in tags]
    times = []
    states, controls = {}, {}
    for j, T in enumerate(trajs):
        ts = np.array(T.getTimeMat()).reshape(-1)
        times.append(ts if j == 0 else ts[1:])
        for n in T.getStateNames():
            v = T.getState(n)
            vals = np.array([v[i] for i in range(v.size())])
            states.setdefault(n, []).append(vals if j == 0 else vals[1:])
        for n in T.getControlNames():
            v = T.getControl(n)
            vals = np.array([v[i] for i in range(v.size())])
            controls.setdefault(n, []).append(vals if j == 0 else vals[1:])
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


def summarize(tags):
    J = 0.0
    for t in tags:
        p = os.path.join(OUT, f"moco_inverse_summary_{t}.txt")
        for line in open(p, encoding="utf-8"):
            if line.startswith("objective="):
                J += float(line.split("=", 1)[1])
                break
    return J


def main():
    global S, GUESS, PREV_J
    MAX_ROUNDS = 30
    PREV_METRICS = np.array([0.890631, 0.741878, 0.402088])
    consecutive_stable = 0
    s_used_last = S
    for r in range(3, MAX_ROUNDS + 1):
        print(f"\n===== outer round {r}: s={S} =====", flush=True)
        tags = []
        with ThreadPoolExecutor(max_workers=4) as ex:
            futs = [ex.submit(solve_chunk, r, ci, S, GUESS) for ci in range(4)]
            for fut in futs:
                tag, rc = fut.result()
                print(f"  {tag}: rc={rc}", flush=True)
                if rc != 0:
                    raise RuntimeError(f"chunk failed: {tag}")
                tags.append(tag)
        tags = [f"B30_right_r{r}_chunk{ci}" for ci in range(4)]
        merged = os.path.join(OUT, f"moco_inverse_solution_B30_right_r{r}_merged.sto")
        T = merge(tags, merged)
        J = summarize(tags)

        ts = np.array(T.getTimeMat()).reshape(-1)
        h = ts[1] - ts[0]
        w = np.ones(len(ts)) * h
        w[0] = w[-1] = h / 2
        s_next, metrics = [], []
        for m, sv in zip(MUSCLES, S):
            u = np.array(T.getControlMat("/forceset/" + m)).reshape(-1)
            e = np.interp(ts, ref_t, ref_col[m])
            s_next.append(float(np.clip(np.sum(w * e * u) / np.sum(w * u * u), 0.5, 2)))
            target = e / sv
            c = np.dot(u, target) / (np.linalg.norm(u) * np.linalg.norm(target))
            rr = np.corrcoef(u, target)[0, 1]
            nr = np.sqrt(np.mean((u - target) ** 2)) / np.sqrt(np.mean(target ** 2))
            metrics.append((c, rr, nr))
        R_s = max(abs(a - b) / b for a, b in zip(s_next, S))
        rJ = abs(J - PREV_J) / PREV_J
        meanm = np.mean(metrics, axis=0)
        metric_change = float(np.max(np.abs(meanm - PREV_METRICS)))
        stable = (R_s < 0.01 and rJ < 0.001 and metric_change < 0.002)
        consecutive_stable = consecutive_stable + 1 if stable else 0

        with open(prog, "a", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow([
                r, f"{J:.6f}", f"{meanm[0]:.4f}", f"{meanm[1]:.4f}",
                f"{meanm[2]:.4f}", *[f"{v:.6f}" for v in S],
                f"{R_s:.6f}", f"{rJ:.6f}",
            ])
        with open(os.path.join(OUT, f"moco_inverse_s_used_B30_right_r{r}.csv"),
                  "w", encoding="utf-8") as f:
            f.write("segment,muscle,s\n")
            for m, v in zip(MUSCLES, S):
                f.write(f"0,{m},{v:.8f}\n")
        print(f"  J={J:.6f} rJ={rJ*100:.3f}% R_s={R_s*100:.2f}% "
              f"mean_cos={meanm[0]:.4f} mean_r={meanm[1]:.4f} "
              f"stable={stable} ({consecutive_stable}/2)", flush=True)
        print(f"  s_next={s_next}", flush=True)
        s_used_last = S
        S = s_next
        GUESS = merged
        PREV_J = J
        PREV_METRICS = meanm
        if consecutive_stable >= 2:
            print(f"  STOP: convergence criteria met at round {r}", flush=True)
            break
    else:
        print(f"  WARNING: stopped at MAX_ROUNDS={MAX_ROUNDS} without 2 stable rounds",
              flush=True)

    # final candidate = last merged
    final = os.path.join(OUT, "moco_inverse_solution_B30_right_final.sto")
    os.replace(GUESS, final)
    with open(os.path.join(OUT, "moco_inverse_s_used_B30_right_final.csv"),
              "w", encoding="utf-8") as f:
        f.write("segment,muscle,s\n")
        for m, v in zip(MUSCLES, s_used_last):
            f.write(f"0,{m},{v:.8f}\n")
    print(f"\nFINAL -> {final}; s_used={s_used_last}; J={PREV_J}", flush=True)


if __name__ == "__main__":
    main()
