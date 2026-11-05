"""Alternating optimization of the reverse-scaled B-form scale factors.

For fixed s, minimize the MocoInverse objective with the reverse-scaled
tracking goal (no free parameters).  Then update each scale factor by the
closed-form block minimizer

    s_i = clip( <e_i, u_i> / <u_i, u_i>,  s_min, s_max )

and repeat with the previous solution as a warm start.  This avoids the
sparsity-detection problem of the joint free-s problem.

Usage:
    python run_alt_s.py --start 199.5 --end 199.8 --mesh 0.02 \
        --weight 20 --tag alt03 --max-iter 6
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
MOCO = os.path.join(HERE, "moco_inverse.py")
REF = os.path.join(HERE, "outputs", "emg_excitation_199p5_204p0.sto")
MUSCLES = [
    "rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r",
    "rect_fem_l", "vasti_l", "hamstrings_l", "tib_ant_l", "gastroc_l",
]


def read_sto(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    i = next(k for k, l in enumerate(lines)
             if l.strip().lower().startswith("endheader"))
    labels = lines[i + 1].strip().split("\t")
    data = np.genfromtxt(path, skip_header=i + 2)
    return labels, data


def update_s(traj_path, ref_path, s_old, lo=0.5, hi=2.0):
    labels, data = read_sto(traj_path)
    t = data[:, 0]
    rlab, rdata = read_sto(ref_path)
    rt = rdata[:, 0]
    s_new = np.array(s_old, float)
    for j, m in enumerate(MUSCLES):
        u = data[:, labels.index(f"/forceset/{m}")]
        e = np.interp(t, rt, rdata[:, rlab.index(f"/forceset/{m}")])
        denom = float(np.dot(u, u))
        if denom <= 0.0:
            s_new[j] = s_old[j]
        else:
            s_new[j] = np.clip(np.dot(e, u) / denom, lo, hi)
    return s_new


def write_guess(src, dst, s_vals):
    # Import lazily: only this helper needs OpenSim.
    sys.path.insert(0, HERE)
    import moco_bootstrap
    moco_bootstrap.ensure_casadi_plugins()
    import opensim as osim
    tr = osim.MocoTrajectory(src)
    for m, v in zip(MUSCLES, s_vals):
        tr.setParameter(f"gain_{m}", float(v))
    tr.write(dst)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=float, default=199.5)
    ap.add_argument("--end", type=float, default=199.8)
    ap.add_argument("--mesh", type=float, default=0.02)
    ap.add_argument("--weight", type=float, default=20.0)
    ap.add_argument("--tag", default="alt")
    ap.add_argument("--max-iter", type=int, default=6)
    ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--s-min", type=float, default=0.5)
    ap.add_argument("--s-max", type=float, default=2.0)
    ap.add_argument("--s-init", type=float, default=1.0)
    args = ap.parse_args(argv)

    out = os.path.join(HERE, "outputs")
    os.makedirs(out, exist_ok=True)
    s = np.full(len(MUSCLES), args.s_init, float)
    prev_sol = None
    for k in range(args.max_iter):
        tag = f"{args.tag}_{k}"
        cmd = [
            sys.executable, MOCO,
            "--start", str(args.start), "--end", str(args.end),
            "--mesh", str(args.mesh),
            "--sagittal-cost", "--nonsagittal-weight", "1e-3",
            "--emg-tracking", REF,
            "--emg-weight", str(args.weight),
            "--emg-rev-goal",
            "--emg-s-values", ",".join(f"{v:.10g}" for v in s),
            "--tag", tag,
        ]
        if prev_sol:
            cmd += ["--emg-initial-guess", prev_sol]
        log = os.path.join(out, f"alt_{args.tag}_{k}.log")
        print(f"\n=== iteration {k}: s = {np.array2string(s, precision=4)}")
        with open(log, "w", encoding="utf-8") as f:
            r = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(f"iteration {k} failed (exit {r.returncode}); see {log}")
            return 1

        seg = os.path.join(out, "segments", "seg_00.sto")
        if not os.path.isfile(seg):
            print(f"missing segment trajectory: {seg}")
            return 1
        s_new = update_s(seg, REF, s, args.s_min, args.s_max)
        dmax = float(np.max(np.abs(s_new - s)))
        print("  s_new = " + ", ".join(
            f"{m}={v:.4f}" for m, v in zip(MUSCLES, s_new)))
        print(f"  max|ds| = {dmax:.6g}")

        # Warm start for the next iteration, with parameters set to s_new.
        guess = os.path.join(out, f"alt_guess_{args.tag}_{k + 1}.sto")
        write_guess(seg, guess, s_new)
        prev_sol = guess
        s = s_new
        if dmax < args.tol:
            print(f"converged after {k + 1} iterations")
            break
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
