"""Long-run B30 right-leg outer loop with checkpoint/resume and progress logs.

This script is intentionally self-contained and portable.  It orchestrates
`experiments/moco/moco_inverse.py` for the four 7.5 s chunks, merges the
retained trajectories, updates the five per-muscle s values, and records
detailed convergence diagnostics.

Two s-update domains are supported:

  --s-domain solve
      Use the 3 s solve windows written by moco_inverse.py (including the
      0.5 s burn-in on each side).  This is consistent with the inner
      objective, which sums the same solve windows.  Recommended.

  --s-domain retained
      Use the stitched 30 s retained trajectory.  This is the domain used by
      the earlier v1 runs; it is NOT the same as the inner objective and must
      be reported as a different (heuristic) update rule.

No optimization is run by --dry-run.
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import opensim as osim

MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r"]
CHUNKS = [(210.0, 217.5), (217.5, 225.0), (225.0, 232.5), (232.5, 240.0)]
S6 = [1.662868, 1.493015, 0.851994, 1.177731, 1.068418]
REF_NAME = "emg_excitation_210_240_right.sto"
DEFAULT_GUESS = "moco_inverse_solution_B30_right_final.sto"


def parse_s(text):
    vals = [float(x) for x in text.replace(",", " ").split()]
    if len(vals) != len(MUSCLES):
        raise argparse.ArgumentTypeError(
            f"need {len(MUSCLES)} s values, got {len(vals)}")
    return vals


def read_ref(outdir):
    path = os.path.join(outdir, REF_NAME)
    lines = open(path, encoding="utf-8").readlines()
    ih = [i for i, line in enumerate(lines)
          if line.strip().lower().startswith("endheader")][0]
    labels = lines[ih + 1].strip().split("\t")
    data = np.genfromtxt(path, skip_header=ih + 2)
    t = data[:, 0]
    cols = {m: data[:, labels.index("/forceset/" + m)] for m in MUSCLES}
    return t, cols


def tracking_obj_T(T, ref_t, ref_cols, s):
    ts = np.array(T.getTimeMat()).reshape(-1)
    h = ts[1] - ts[0]
    w = np.ones(len(ts)) * h
    w[0] = w[-1] = h / 2
    total = 0.0
    for m, sv in zip(MUSCLES, s):
        u = np.array(T.getControlMat("/forceset/" + m)).reshape(-1)
        e = np.interp(ts, ref_t, ref_cols[m])
        total += 20.0 * np.sum(w * (e - sv * u) ** 2)
    return float(total)


def tracking_obj_files(files, ref_t, ref_cols, s):
    total = 0.0
    for path in files:
        T = osim.MocoTrajectory(path)
        total += tracking_obj_T(T, ref_t, ref_cols, s)
    return float(total)


def metric_for(T, ref_t, ref_cols, s):
    ts = np.array(T.getTimeMat()).reshape(-1)
    h = ts[1] - ts[0]
    w = np.ones(len(ts)) * h
    w[0] = w[-1] = h / 2
    out = []
    for m, sv in zip(MUSCLES, s):
        u = np.array(T.getControlMat("/forceset/" + m)).reshape(-1)
        e = np.interp(ts, ref_t, ref_cols[m])
        target = e / sv
        c = np.dot(u, target) / (np.linalg.norm(u) * np.linalg.norm(target))
        r = np.corrcoef(u, target)[0, 1]
        nr = np.sqrt(np.mean((u - target) ** 2)) / np.sqrt(np.mean(target ** 2))
        out.append((c, r, nr))
    return np.asarray(out)


def s_hat_retained(T, ref_t, ref_cols):
    ts = np.array(T.getTimeMat()).reshape(-1)
    h = ts[1] - ts[0]
    w = np.ones(len(ts)) * h
    w[0] = w[-1] = h / 2
    s = []
    for m in MUSCLES:
        u = np.array(T.getControlMat("/forceset/" + m)).reshape(-1)
        e = np.interp(ts, ref_t, ref_cols[m])
        s.append(float(np.clip(np.sum(w * e * u) / np.sum(w * u * u), 0.5, 2)))
    return s


def s_hat_solve(seg_files, ref_t, ref_cols):
    """Integrate over every 3 s solve window, matching the inner objective."""
    num = np.zeros(len(MUSCLES))
    den = np.zeros(len(MUSCLES))
    for path in seg_files:
        T = osim.MocoTrajectory(path)
        ts = np.array(T.getTimeMat()).reshape(-1)
        h = ts[1] - ts[0]
        w = np.ones(len(ts)) * h
        w[0] = w[-1] = h / 2
        for k, m in enumerate(MUSCLES):
            u = np.array(T.getControlMat("/forceset/" + m)).reshape(-1)
            e = np.interp(ts, ref_t, ref_cols[m])
            num[k] += np.sum(w * e * u)
            den[k] += np.sum(w * u * u)
    return [float(np.clip(a / b, 0.5, 2)) for a, b in zip(num, den)]


def merge(tags, outdir, out_path):
    trajs = [osim.MocoTrajectory(
        os.path.join(outdir, f"moco_inverse_solution_{t}.sto")) for t in tags]
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


def parse_components(log_path):
    """Return (effort, tracking, number of breakdown blocks) from a chunk log."""
    txt = open(log_path, encoding="utf-8", errors="ignore").read()
    eff = trk = 0.0
    n = 0
    pat = re.compile(
        r"Breakdown of objective \(including weights\):\s*\n"
        r"\s*excitation_effort:\s*([0-9.eE+-]+)\s*\n"
        r"\s*emg_tracking:\s*([0-9.eE+-]+)")
    for m in pat.finditer(txt):
        eff += float(m.group(1))
        trk += float(m.group(2))
        n += 1
    return eff, trk, n


def max_activation_change(prev, new):
    vals = []
    for n in prev.getStateNames():
        if not n.endswith("/activation"):
            continue
        a = np.array([prev.getState(n)[i] for i in range(prev.getState(n).size())])
        b = np.array([new.getState(n)[i] for i in range(new.getState(n).size())])
        vals.append(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-8))
    return float(max(vals)) if vals else float("nan")


def max_u_change(prev, new):
    vals = []
    for m in MUSCLES:
        a = np.array(prev.getControlMat("/forceset/" + m)).reshape(-1)
        b = np.array(new.getControlMat("/forceset/" + m)).reshape(-1)
        vals.append(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-8))
    return float(max(vals))


def saturation(new):
    out = {}
    for n in new.getStateNames():
        if not n.endswith("/activation"):
            continue
        v = np.array([new.getState(n)[i] for i in range(new.getState(n).size())])
        out[n.split("/")[2]] = (float(v.max()), float(np.mean(v > 0.99)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--init-s", type=parse_s, default=S6,
                    help="initial five s values (default: r7 s_used)")
    ap.add_argument("--init-guess", default=DEFAULT_GUESS,
                    help="initial 30 s MocoTrajectory used as guess")
    ap.add_argument("--s-domain", choices=["solve", "retained"], default="solve",
                    help="s-update domain; solve matches the inner objective")
    ap.add_argument("--max-updates", type=int, default=10,
                    help="updates in this invocation")
    ap.add_argument("--max-total-updates", type=int, default=40,
                    help="stop after this many total effective updates")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--mesh", type=float, default=0.02)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--max-iter", type=int, default=5000)
    ap.add_argument("--resume", default=None,
                    help="checkpoint.json to resume from")
    ap.add_argument("--outdir", default=None,
                    help="experiments/moco/outputs by default")
    ap.add_argument("--tag-prefix", default="B30_cont")
    ap.add_argument("--dry-run", action="store_true",
                    help="print commands and exit without solving")
    args = ap.parse_args()
    if args.workers > len(CHUNKS):
        print(f"NOTE: only {len(CHUNKS)} chunks exist; capping --workers "
              f"from {args.workers} to {len(CHUNKS)}", flush=True)
        args.workers = len(CHUNKS)

    root = os.path.dirname(os.path.abspath(__file__))
    outdir = args.outdir or os.path.join(root, "experiments", "moco", "outputs")
    logdir = os.path.join(outdir, "continuous_logs")
    os.makedirs(logdir, exist_ok=True)
    py = sys.executable

    s_current = list(args.init_s)
    guess = os.path.abspath(os.path.join(outdir, args.init_guess))
    update_id = 0
    total_updates = 0
    resume_stable = 0
    resume_J = None

    ref_path = os.path.join(outdir, REF_NAME)
    if not os.path.isfile(ref_path):
        print(f"ERROR: reference not found: {ref_path}", flush=True)
        print("Copy emg_excitation_210_240_right.sto into outputs/ first.", flush=True)
        return 2

    if args.resume and not os.path.isfile(args.resume):
        print(f"ERROR: checkpoint not found: {args.resume}", flush=True)
        print("A checkpoint is created only after the first successful update.", flush=True)
        print("Start WITHOUT --resume for the first run, e.g.:", flush=True)
        print("  python .\\run_B30_continue.py --s-domain solve "
              "--max-updates 1 --max-total-updates 40 --workers 4", flush=True)
        return 2

    if not args.resume and not os.path.isfile(guess):
        print(f"ERROR: initial guess not found: {guess}", flush=True)
        print("Copy moco_inverse_solution_B30_right_final.sto into outputs/ first.", flush=True)
        return 2

    ref_t, ref_cols = read_ref(outdir)

    if args.resume:
        try:
            ck = json.load(open(args.resume, encoding="utf-8"))
        except Exception as exc:
            print(f"ERROR: cannot read checkpoint {args.resume}: {exc}", flush=True)
            return 2
        s_current = [float(x) for x in ck["s_current"]]
        guess = ck["solution"]
        update_id = int(ck["update_id"])
        total_updates = int(ck["total_updates"])
        pending_s = [float(x) for x in ck["s_hat"]]
        resume_stable = int(ck.get("stable_count", 0))
        resume_J = float(ck["J"]) if "J" in ck else None
        print(f"RESUME {args.resume}: update_id={update_id} s={s_current}", flush=True)
    else:
        pending_s = None

    ckpt_path = os.path.join(outdir, f"{args.tag_prefix}_checkpoint.json")
    progress_path = os.path.join(outdir, f"{args.tag_prefix}_progress.csv")
    if not os.path.isfile(progress_path):
        with open(progress_path, "w", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow([
                "update_id", "timestamp", "J", "effort", "tracking",
                "mean_cos", "mean_r", "mean_nRMSE", "s0", "s1", "s2", "s3", "s4",
                "s_hat0", "s_hat1", "s_hat2", "s_hat3", "s_hat4",
                "R_s", "rJ", "G_s", "u_rel", "act_rel", "stable_count", "bound",
                "solver_ok", "elapsed_s",
            ])

    print("=" * 78, flush=True)
    print(f"B30 continuation | domain={args.s_domain} | workers={args.workers}", flush=True)
    print(f"outdir={outdir}", flush=True)
    print(f"start s={s_current}", flush=True)
    print(f"guess={guess}", flush=True)
    print("=" * 78, flush=True)

    stable_count = resume_stable
    prev_J = resume_J
    prev_metrics = None
    s_last = list(s_current)
    failed = False
    for local in range(args.max_updates):
        if total_updates >= args.max_total_updates:
            print("STOP: max_total_updates reached", flush=True)
            break
        update_id += 1
        total_updates += 1
        t0 = time.time()
        print(f"\n--- update {update_id} / total {total_updates} ---", flush=True)
        print(f"s_prev (last solved) = {['%.6f' % x for x in s_current]}", flush=True)

        # Compute the candidate s BEFORE solving.  On resume use the stored
        # candidate so that a restart does not repeat a parameter update.
        T_cur = osim.MocoTrajectory(guess)
        if pending_s is not None:
            s_hat = pending_s
            pending_s = None
            print("s_hat_from_prev (candidate from checkpoint)", flush=True)
        elif args.s_domain == "retained":
            s_hat = s_hat_retained(T_cur, ref_t, ref_cols)
        else:
            # solve-domain requires per-segment files from the previous solve.
            # For the first invocation (or after deleting segments), use the
            # retained estimate as a clearly labelled bootstrap candidate.
            s_hat = s_hat_retained(T_cur, ref_t, ref_cols)
            print("s_hat bootstrap=retained (first solve-domain update)", flush=True)
        print(f"s_hat  = {['%.6f' % x for x in s_hat]}", flush=True)

        s_solve = s_hat
        s_last = s_solve
        print(f"s_used (this update) = {['%.6f' % x for x in s_solve]}", flush=True)
        tags = []
        for ci in range(4):
            tags.append(f"{args.tag_prefix}_u{update_id:03d}_c{ci}")

        commands = []
        for ci, (start, end) in enumerate(CHUNKS):
            tag = tags[ci]
            segdir = os.path.join(outdir, f"segments_{tag}")
            cmd = [
                py, os.path.join(root, "experiments", "moco", "moco_inverse.py"),
                "--start", str(start), "--end", str(end),
                "--segment", "2.0", "--overlap", "0.5", "--mesh", str(args.mesh),
                "--sagittal-cost", "--nonsagittal-weight", "1e-3",
                "--emg-tracking", os.path.join(outdir, REF_NAME),
                "--emg-weight", "20", "--emg-leg", "right",
                "--emg-native-b-fixed",
                "--emg-s-values", ",".join(f"{v:.6f}" for v in s_solve),
                "--emg-initial-guess", guess,
                "--emg-sparsity-random", "--tol", str(args.tol),
                "--max-iter", str(args.max_iter),
                "--seg-dir", segdir,
                "--tag", tag,
            ]
            commands.append((cmd, os.path.join(logdir, tag + ".log"), tag))

        if args.dry_run:
            for cmd, log, _ in commands:
                print("DRY:", " ".join(cmd), flush=True)
            continue

        def run_one(item):
            cmd, log, tag = item
            with open(log, "w", encoding="utf-8") as fh:
                p = subprocess.run(cmd, cwd=root, stdout=fh,
                                   stderr=subprocess.STDOUT)
            return tag, p.returncode, log

        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            results = list(ex.map(run_one, commands))
        ok = all(rc == 0 for _, rc, _ in results)
        print("chunks: " + ", ".join(f"{t}:rc={rc}" for t, rc, _ in results), flush=True)
        if not ok:
            print("ERROR: at least one chunk failed; checkpoint not advanced", flush=True)
            for tag, rc, log in results:
                if rc == 0:
                    continue
                print(f"--- tail {os.path.basename(log)} (rc={rc}) ---", flush=True)
                try:
                    lines = open(log, encoding="utf-8", errors="ignore").read().splitlines()
                    for line in lines[-40:]:
                        print("  " + line, flush=True)
                except Exception as exc:
                    print(f"  cannot read log: {exc}", flush=True)
            failed = True
            break

        merged = os.path.join(outdir, f"moco_inverse_solution_{args.tag_prefix}_u{update_id:03d}_merged.sto")
        T_new = merge(tags, outdir, merged)
        metrics = metric_for(T_new, ref_t, ref_cols, s_solve)
        meanm = metrics.mean(axis=0)
        J = 0.0
        effort = tracking = 0.0
        nblocks = 0
        for _, log, _ in commands:
            e, tr, n = parse_components(log)
            effort += e
            tracking += tr
            nblocks += n
            # summary objective for this chunk
            tag = os.path.basename(log)[:-4]
            sp = os.path.join(outdir, f"moco_inverse_summary_{tag}.txt")
            for line in open(sp, encoding="utf-8"):
                if line.startswith("objective="):
                    J += float(line.split("=", 1)[1])
                    break
        if args.s_domain == "solve":
            seg_files = []
            for ci in range(4):
                d = os.path.join(outdir, f"segments_{tags[ci]}")
                seg_files += sorted(
                    os.path.join(d, f) for f in os.listdir(d)
                    if f.startswith("seg_") and f.endswith(".sto"))
            s_hat_new = s_hat_solve(seg_files, ref_t, ref_cols)
            G_s = (tracking_obj_files(seg_files, ref_t, ref_cols, s_solve)
                   - tracking_obj_files(seg_files, ref_t, ref_cols, s_hat_new))
        else:
            s_hat_new = s_hat_retained(T_new, ref_t, ref_cols)
            G_s = (tracking_obj_T(T_new, ref_t, ref_cols, s_solve)
                   - tracking_obj_T(T_new, ref_t, ref_cols, s_hat_new))
        R_s = max(abs(a - b) / b for a, b in zip(s_hat_new, s_solve))
        rJ = (abs(J - prev_J) / prev_J) if prev_J else float("nan")
        u_rel = max_u_change(T_cur, T_new)
        act_rel = max_activation_change(T_cur, T_new)
        bound = [m for m, v in zip(MUSCLES, s_solve)
                 if v <= 0.5005 or v >= 1.9995]
        sat = saturation(T_new)
        stable = (R_s < 0.01 and u_rel < 0.01 and act_rel < 0.02
                  and (np.isnan(rJ) or rJ < 0.001))
        stable_count = stable_count + 1 if stable else 0
        elapsed = time.time() - t0
        print(f"J={J:.6f} effort={effort:.6f} tracking={tracking:.6f} "
              f"blocks={nblocks}", flush=True)
        print(f"mean cos={meanm[0]:.4f} r={meanm[1]:.4f} nRMSE={meanm[2]:.4f}", flush=True)
        print(f"R_s={R_s*100:.3f}% rJ={rJ*100 if not np.isnan(rJ) else float('nan'):.3f}% "
              f"G_s={G_s:.6f} u_rel={u_rel*100:.2f}% act_rel={act_rel*100:.2f}%", flush=True)
        print(f"stable={stable} ({stable_count}/3) bound={bound} elapsed={elapsed:.1f}s",
              flush=True)

        with open(progress_path, "a", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow([
                update_id, time.strftime("%Y-%m-%dT%H:%M:%S"),
                f"{J:.6f}", f"{effort:.6f}", f"{tracking:.6f}",
                f"{meanm[0]:.4f}", f"{meanm[1]:.4f}", f"{meanm[2]:.4f}",
                *[f"{v:.6f}" for v in s_solve],
                *[f"{v:.6f}" for v in s_hat_new],
                f"{R_s:.6f}", f"{rJ:.6f}", f"{G_s:.6f}", f"{u_rel:.6f}", f"{act_rel:.6f}",
                stable_count, ";".join(bound), ok, f"{elapsed:.1f}",
            ])
        ck = {
            "schema": 1,
            "domain": args.s_domain,
            "update_id": update_id,
            "total_updates": total_updates,
            "s_current": s_solve,
            "s_hat": s_hat_new,
            "solution": merged,
            "guess": guess,
            "J": J,
            "effort": effort,
            "tracking": tracking,
            "R_s": R_s,
            "rJ": rJ,
            "G_s": G_s,
            "u_rel": u_rel,
            "act_rel": act_rel,
            "stable_count": stable_count,
            "bound": bound,
            "saturation": sat,
            "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        }
        tmp = ckpt_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(ck, f, indent=2)
        os.replace(tmp, ckpt_path)
        s_current = s_hat_new
        guess = merged
        pending_s = s_hat_new
        prev_J = J
        prev_metrics = meanm
        if stable_count >= 3:
            print("STOP: 3 consecutive stable updates", flush=True)
            break

    final = os.path.join(outdir, f"{args.tag_prefix}_final.sto")
    if not args.dry_run and not failed and os.path.isfile(guess):
        shutil.copyfile(guess, final)
        with open(os.path.join(outdir, f"{args.tag_prefix}_final_s_used.csv"),
                  "w", encoding="utf-8") as f:
            f.write("segment,muscle,s\n")
            for m, v in zip(MUSCLES, s_last):
                f.write(f"0,{m},{v:.8f}\n")
        print(f"\nFINAL -> {final}", flush=True)
        print(f"final s_used -> {os.path.join(outdir, f'{args.tag_prefix}_final_s_used.csv')}", flush=True)
        print(f"checkpoint -> {ckpt_path}", flush=True)
        print(f"progress   -> {progress_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
