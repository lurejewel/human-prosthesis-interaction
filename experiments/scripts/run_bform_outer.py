"""B-form inner/outer alternating MocoInverse driver for one JW trial folder.

Generation 0 (bootstrap) solves with fixed s=1 for all 10 measured muscles and
no initial guess.  Every later generation fixes the closed-form s update from
the previous trajectory, uses that trajectory as the initial guess, and solves
independent chunks in parallel.  The outer update is the least-squares block
minimizer s = sum(w*e*u)/sum(w*u*u), clipped to [0.5, 2].

Typical use from a trial folder:
    python run_bform_outer.py --workers 4 --bootstrap-only
    python run_bform_outer.py --resume bform_outputs/checkpoint.json --workers 4

Or run bootstrap and all outer updates automatically:
    python run_bform_outer.py --workers 4 --max-updates 40

No optimization is run by --dry-run.
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moco_bootstrap
moco_bootstrap.ensure_casadi_plugins()
import opensim as osim

MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r",
           "rect_fem_l", "vasti_l", "hamstrings_l", "tib_ant_l", "gastroc_l"]
S_MIN, S_MAX = 0.5, 2.0


def one(folder, pattern):
    import glob
    hits = sorted(glob.glob(os.path.join(folder, pattern)))
    if not hits:
        raise SystemExit(f"no file matching {pattern} in {folder}")
    return hits[0]


def read_ref(path):
    lines = open(path, encoding="utf-8").readlines()
    ih = next(i for i, line in enumerate(lines)
              if line.strip().lower().startswith("endheader"))
    labels = lines[ih + 1].strip().split("\t")
    data = np.genfromtxt(path, skip_header=ih + 2)
    t = data[:, 0]
    cols = {m: data[:, labels.index("/forceset/" + m)] for m in MUSCLES}
    return t, cols


def trajectory_values(T, name, control=True):
    if control:
        return np.array(T.getControlMat(name)).reshape(-1)
    v = T.getState(name)
    return np.array([v[i] for i in range(v.size())])


def s_hat_from_files(files, ref_t, ref_cols, s_old):
    num = np.zeros(len(MUSCLES))
    den = np.zeros(len(MUSCLES))
    for path in files:
        T = osim.MocoTrajectory(path)
        ts = np.array(T.getTimeMat()).reshape(-1)
        h = float(np.median(np.diff(ts)))
        w = np.full(len(ts), h)
        w[0] = w[-1] = h / 2
        for k, m in enumerate(MUSCLES):
            u = trajectory_values(T, "/forceset/" + m)
            e = np.interp(ts, ref_t, ref_cols[m])
            num[k] += np.sum(w * e * u)
            den[k] += np.sum(w * u * u)
    return [float(np.clip(num[k] / den[k], S_MIN, S_MAX)) if den[k] > 1e-12
            else float(s_old[k]) for k in range(len(MUSCLES))]


def merge_trajectories(paths, out_path):
    trajs = [osim.MocoTrajectory(path) for path in paths]
    times, states, controls = [], {}, {}
    for j, T in enumerate(trajs):
        ts = np.array(T.getTimeMat()).reshape(-1)
        times.append(ts if j == 0 else ts[1:])
        for name in T.getStateNames():
            vals = trajectory_values(T, name, control=False)
            states.setdefault(name, []).append(vals if j == 0 else vals[1:])
        for name in T.getControlNames():
            vals = trajectory_values(T, name)
            controls.setdefault(name, []).append(vals if j == 0 else vals[1:])
    t = np.concatenate(times)
    merged = trajs[0]
    merged.setNumTimes(len(t))
    merged.setTime(t.tolist())
    for name in merged.getStateNames():
        merged.setState(name, np.concatenate(states[name]).tolist())
    for name in merged.getControlNames():
        merged.setControl(name, np.concatenate(controls[name]).tolist())
    merged.write(out_path)
    return merged


def max_u_change(prev, new):
    vals = []
    for m in MUSCLES:
        a = trajectory_values(prev, "/forceset/" + m)
        b = trajectory_values(new, "/forceset/" + m)
        vals.append(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-8))
    return float(max(vals))


def max_activation_change(prev, new):
    vals = []
    for name in prev.getStateNames():
        if not name.endswith("/activation"):
            continue
        a = trajectory_values(prev, name, control=False)
        b = trajectory_values(new, name, control=False)
        vals.append(np.linalg.norm(b - a) / max(np.linalg.norm(a), 1e-8))
    return float(max(vals)) if vals else float("nan")


def objective_from_summaries(paths):
    total = 0.0
    for path in paths:
        for line in open(path, encoding="utf-8"):
            if line.startswith("objective="):
                total += float(line.split("=", 1)[1])
                break
    return total


def parse_s(text):
    if text is None:
        return [1.0] * len(MUSCLES)
    vals = [float(x) for x in text.replace(",", " ").split()]
    if len(vals) != len(MUSCLES):
        raise argparse.ArgumentTypeError(
            f"need {len(MUSCLES)} s values, got {len(vals)}")
    return vals


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunks", type=int, default=4)
    ap.add_argument("--mesh", type=float, default=0.002)
    ap.add_argument("--segment", type=float, default=2.0)
    ap.add_argument("--overlap", type=float, default=0.5)
    ap.add_argument("--emg-weight", type=float, default=20.0)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--max-iter", type=int, default=5000)
    ap.add_argument("--max-updates", type=int, default=40,
                    help="maximum outer updates after bootstrap")
    ap.add_argument("--stable-updates", type=int, default=3)
    ap.add_argument("--s-domain", choices=["solve", "retained"], default="solve")
    ap.add_argument("--init-s", type=parse_s, default=None,
                    help="10 fixed s values for bootstrap (default: all 1)")
    ap.add_argument("--init-guess", default=None,
                    help="optional B-form trajectory used only for the bootstrap solve")
    ap.add_argument("--bootstrap-only", action="store_true",
                    help="stop after generation 0; resume later from checkpoint.json")
    ap.add_argument("--resume", default=None,
                    help="bform_outputs/checkpoint.json from a previous run")
    ap.add_argument("--tag-prefix", default="BFORM")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    folder = os.path.abspath(args.data_dir)
    os.chdir(folder)
    summary = json.load(open(os.path.join(folder, "steady_state_summary.json")))
    a0, a1 = summary["analysis_window"]
    edges = np.linspace(a0, a1, args.chunks + 1)
    chunks = [(float(edges[i]), float(edges[i + 1])) for i in range(args.chunks)]
    workers = min(args.workers, args.chunks)
    ref_path = one(folder, "*_emg_reference.sto")
    model_path = one(folder, "*_simbody.osim")
    ik_path = one(folder, "*_ik_filtered.mot")
    grf_path = one(folder, "*_grf_filtered.xml")
    moco_path = os.path.join(folder, "moco_inverse.py")
    if not os.path.isfile(moco_path):
        raise SystemExit(f"missing {moco_path}")
    outdir = os.path.join(folder, "bform_outputs")
    ckpt_path = os.path.join(outdir, "checkpoint.json")
    progress_path = os.path.join(outdir, "progress.csv")

    update_id = 0
    s_solve = list(args.init_s) if args.init_s is not None else [1.0] * len(MUSCLES)
    guess = os.path.abspath(args.init_guess) if args.init_guess else None
    prev_J = None
    stable_count = 0
    if args.resume:
        if not os.path.isfile(args.resume):
            raise SystemExit(f"checkpoint not found: {args.resume}")
        ck = json.load(open(args.resume, encoding="utf-8"))
        update_id = int(ck["next_update"])
        s_solve = [float(v) for v in ck["s_next"]]
        guess = os.path.abspath(ck["solution"])
        prev_J = float(ck["J"])
        stable_count = int(ck.get("stable_count", 0))
    if guess and not os.path.isfile(guess):
        raise SystemExit(f"initial guess not found: {guess}")

    print(f"trial      : {os.path.basename(folder)}", flush=True)
    print(f"window     : {a0:.3f}-{a1:.3f} s; chunks={chunks}", flush=True)
    print(f"B-form     : 10 muscles, weight={args.emg_weight}, mesh={args.mesh}", flush=True)
    print(f"start      : update={update_id}, s={np.round(s_solve, 6).tolist()}", flush=True)
    print(f"guess      : {guess or 'none (bootstrap)'}", flush=True)
    print(f"parallel   : {workers} workers; domain={args.s_domain}", flush=True)

    def commands_for(current_s, current_guess, current_update):
        commands = []
        for ci, (start, end) in enumerate(chunks):
            tag = f"{args.tag_prefix}_u{current_update:03d}_c{ci}"
            chunk_out = os.path.join(outdir, "chunks", tag)
            seg_dir = os.path.join(chunk_out, "segments")
            cmd = [
                sys.executable, moco_path,
                "--start", f"{start:.6f}", "--end", f"{end:.6f}",
                "--model-file", model_path, "--ik-file", ik_path,
                "--grf-xml", grf_path, "--out-dir", chunk_out,
                "--segment", str(args.segment), "--overlap", str(args.overlap),
                "--mesh", str(args.mesh),
                "--sagittal-cost", "--nonsagittal-weight", "1e-3",
                "--emg-tracking", ref_path,
                "--emg-weight", str(args.emg_weight), "--emg-leg", "both",
                "--emg-native-b-fixed",
                "--emg-s-values", ",".join(f"{v:.8f}" for v in current_s),
                "--emg-sparsity-random", "--tol", str(args.tol),
                "--max-iter", str(args.max_iter),
                "--seg-dir", seg_dir, "--tag", tag,
            ]
            if current_guess:
                cmd += ["--emg-initial-guess", current_guess]
            commands.append((tag, chunk_out, cmd))
        return commands

    if args.dry_run:
        for _, _, cmd in commands_for(s_solve, guess, update_id):
            print("DRY:", " ".join(cmd), flush=True)
        return 0

    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    if not os.path.isfile(progress_path):
        with open(progress_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["update", "timestamp", "J", "R_s", "rJ", "u_rel",
                             "act_rel", "stable_count", "bound", "elapsed_s",
                             *[f"s_used_{m}" for m in MUSCLES],
                             *[f"s_next_{m}" for m in MUSCLES]])
    ref_t, ref_cols = read_ref(ref_path)

    while True:
        if update_id > args.max_updates:
            print("STOP: maximum outer updates reached", flush=True)
            break
        t0 = time.time()
        print(f"\n--- update {update_id} | s_used={np.round(s_solve, 6).tolist()} ---",
              flush=True)
        commands = commands_for(s_solve, guess, update_id)

        def run_one(item):
            tag, chunk_out, cmd = item
            os.makedirs(chunk_out, exist_ok=True)
            log_path = os.path.join(outdir, "logs", tag + ".log")
            with open(log_path, "w", encoding="utf-8") as fh:
                proc = subprocess.run(cmd, cwd=folder, stdout=fh,
                                      stderr=subprocess.STDOUT)
            return tag, proc.returncode, log_path

        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(run_one, commands))
        print("chunks: " + ", ".join(f"{tag}:rc={rc}" for tag, rc, _ in results),
              flush=True)
        if any(rc != 0 for _, rc, _ in results):
            print("ERROR: chunk failed; checkpoint not advanced", flush=True)
            return 1

        solution_paths = [os.path.join(chunk_out, f"moco_inverse_solution_{tag}.sto")
                          for tag, chunk_out, _ in commands]
        merged_path = os.path.join(
            outdir, f"{args.tag_prefix}_u{update_id:03d}_merged.sto")
        T_new = merge_trajectories(solution_paths, merged_path)
        if args.s_domain == "solve":
            update_files = []
            for _, chunk_out, _ in commands:
                seg_dir = os.path.join(chunk_out, "segments")
                update_files += sorted(
                    os.path.join(seg_dir, f) for f in os.listdir(seg_dir)
                    if f.startswith("seg_") and f.endswith(".sto"))
        else:
            update_files = [merged_path]
        s_next = s_hat_from_files(update_files, ref_t, ref_cols, s_solve)
        R_s = max(abs(a - b) / b for a, b in zip(s_next, s_solve))
        J = objective_from_summaries([
            os.path.join(chunk_out, f"moco_inverse_summary_{tag}.txt")
            for tag, chunk_out, _ in commands])
        rJ = abs(J - prev_J) / prev_J if prev_J else float("nan")
        if guess:
            T_prev = osim.MocoTrajectory(guess)
            u_rel = max_u_change(T_prev, T_new)
            act_rel = max_activation_change(T_prev, T_new)
        else:
            u_rel = act_rel = float("nan")
        stable = (update_id > 0 and R_s < 0.01 and u_rel < 0.01
                  and act_rel < 0.02 and (np.isnan(rJ) or rJ < 0.001))
        stable_count = stable_count + 1 if stable else 0
        bound = [m for m, v in zip(MUSCLES, s_solve)
                 if v <= S_MIN + 0.0005 or v >= S_MAX - 0.0005]
        elapsed = time.time() - t0
        print(f"J={J:.6f} R_s={R_s:.4%} rJ={rJ:.4%} "
              f"u_rel={u_rel:.4%} act_rel={act_rel:.4%}", flush=True)
        print(f"s_next={np.round(s_next, 6).tolist()}", flush=True)
        print(f"stable={stable} ({stable_count}/{args.stable_updates}); "
              f"bound={bound}; elapsed={elapsed:.1f}s", flush=True)

        with open(progress_path, "a", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow([
                update_id, time.strftime("%Y-%m-%dT%H:%M:%S"), f"{J:.6f}",
                f"{R_s:.6f}", f"{rJ:.6f}", f"{u_rel:.6f}", f"{act_rel:.6f}",
                stable_count, ";".join(bound), f"{elapsed:.1f}",
                *[f"{v:.8f}" for v in s_solve],
                *[f"{v:.8f}" for v in s_next]])
        ck = {"next_update": update_id + 1, "s_next": s_next,
              "s_used": s_solve, "solution": merged_path, "J": J,
              "stable_count": stable_count, "domain": args.s_domain,
              "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S")}
        tmp = ckpt_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(ck, fh, indent=2)
        os.replace(tmp, ckpt_path)

        guess = merged_path
        s_solve = s_next
        prev_J = J
        update_id += 1
        if args.bootstrap_only:
            print(f"\nBootstrap complete -> {merged_path}", flush=True)
            print(f"Resume with: python run_bform_outer.py --resume {ckpt_path} "
                  f"--workers {workers}", flush=True)
            return 0
        if stable_count >= args.stable_updates:
            print("STOP: consecutive stable updates reached", flush=True)
            break
        if update_id > args.max_updates:
            print("STOP: maximum outer updates reached", flush=True)
            break

    final_path = os.path.join(outdir, f"{args.tag_prefix}_final.sto")
    import shutil
    shutil.copyfile(guess, final_path)
    final_s = os.path.join(outdir, f"{args.tag_prefix}_final_s_used.csv")
    with open(final_s, "w", encoding="utf-8") as fh:
        fh.write("muscle,s\n")
        for m, v in zip(MUSCLES, json.load(open(ckpt_path))["s_used"]):
            fh.write(f"{m},{v:.8f}\n")
    print(f"\nFINAL -> {final_path}\nfinal s -> {final_s}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())