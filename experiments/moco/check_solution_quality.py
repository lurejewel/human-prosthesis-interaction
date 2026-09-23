"""Solution-quality check for MocoInverse results.

Reports, for a solution file:

* **mesh-scale chatter** of the muscle excitations and activations.  This is
  the key check.  MocoInverse's cost penalises only the *magnitude* of the
  excitation, and the collocation of the activation dynamics admits a
  near-alternating (zigzag) mode when the mesh interval h is not small
  compared with the activation time constant tau (~10 ms for
  DeGrooteFregly2016).  At h = 20 ms (h/tau = 2) the activation itself picks up
  several percent of ripple, which is not physiological and makes the muscle
  forces unreliable; at h = 2 ms it is essentially gone.

  Metric:  alt_i = y_i - (y_{i-1} + y_{i+1})/2  on the *collocation mesh*
  (the solution is written at mesh/2 spacing, so the mesh grid is every 2nd
  sample).  Reported as RMS(alt)/std(y):  ~0 for a smooth trajectory,
  O(1) when the trajectory alternates every mesh point.

* actuator saturation (muscles pinned at activation 1, residuals pinned at
  their control bound).

Usage
-----
    python check_solution_quality.py                       # default solution
    python check_solution_quality.py outputs/moco_inverse_solution_fine.sto
    python check_solution_quality.py FILE --mesh 0.002
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moco_bootstrap                                     # noqa: E402
moco_bootstrap.ensure_casadi_plugins()
import opensim as osim                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Rough guidance derived from measurements on this dataset:
ALT_WARN = 0.05      # RMS(alt)/std above which a trajectory is "wiggly"
ALT_BAD = 0.15


def alt_ratio(y):
    """RMS of the local alternating component / std of the signal."""
    if len(y) < 5:
        return float("nan")
    alt = y[1:-1] - 0.5 * (y[:-2] + y[2:])
    s = float(np.std(y))
    return float(np.sqrt(np.mean(alt ** 2)) / s) if s > 0 else float("nan")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("solution", nargs="?",
                   default=os.path.join(HERE, "outputs", "moco_inverse_solution.sto"))
    p.add_argument("--mesh", type=float, default=None,
                   help="collocation mesh interval [s]; default: infer as "
                        "2 x the solution's output spacing")
    p.add_argument("--window", type=float, nargs=2, default=None,
                   metavar=("T0", "T1"), help="restrict analysis to a time window")
    args = p.parse_args(argv)

    if not os.path.exists(args.solution):
        sys.exit(f"not found: {args.solution}")

    sol = osim.MocoTrajectory(args.solution)
    n = sol.getNumTimes()
    t = np.array([sol.getTime()[i] for i in range(n)])
    st, ct = sol.exportToStatesTable(), sol.exportToControlsTable()
    dt = float(np.median(np.diff(t)))
    mesh = args.mesh if args.mesh else 2.0 * dt
    stride = max(int(round(mesh / dt)), 1)

    sel = np.ones(n, dtype=bool)
    if args.window:
        sel = (t >= args.window[0]) & (t <= args.window[1])
    idx = np.where(sel)[0][::stride]

    print(f"file        : {args.solution}")
    print(f"output grid : dt = {dt:.4f} s, {n} points, "
          f"{t[0]:.2f}..{t[-1]:.2f} s")
    print(f"mesh        : {mesh:.4f} s (stride {stride}, {len(idx)} mesh samples)")

    muscles, res_ctrl = [], []
    for lab in ct.getColumnLabels():
        short = lab.rsplit("/", 1)[-1]
        if short.startswith(("residual_", "reserve_")):
            res_ctrl.append(lab)
        else:
            muscles.append(short)

    print(f"\n{'muscle':<14}{'alt/std x':>11}{'alt/std a':>11}"
          f"{'alt RMS a':>11}{'peak a':>9}  verdict")
    wx, wa, worst = [], [], []
    for m in muscles:
        a = st.getDependentColumn(f"/forceset/{m}/activation").to_numpy()[idx]
        x = ct.getDependentColumn(f"/forceset/{m}").to_numpy()[idx]
        axr, aar = alt_ratio(x), alt_ratio(a)
        aalt = float(np.sqrt(np.mean((a[1:-1] - 0.5 * (a[:-2] + a[2:])) ** 2)))
        wx.append(axr)
        wa.append(aar)
        worst.append((aar, m))
        verdict = ("OK" if aar < ALT_WARN else
                   "wiggly" if aar < ALT_BAD else "CHATTERING")
        print(f"{m:<14}{axr:>11.3f}{aar:>11.3f}{aalt:>11.4f}"
              f"{a.max():>9.3f}  {verdict}")

    print(f"\n{'MEAN':<14}{np.mean(wx):>11.3f}{np.mean(wa):>11.3f}")
    worst.sort(reverse=True)
    print(f"worst activation chatter: {worst[0][1]} ({worst[0][0]:.3f})")

    print(f"\nactivation at bound (>0.98): "
          f"{sum(1 for m in muscles if st.getDependentColumn(f'/forceset/{m}/activation').to_numpy().max() > 0.98)}"
          f" / {len(muscles)} muscles")
    n_near = 0
    for lab in res_ctrl:
        if lab.rsplit('/', 1)[-1].startswith("residual_"):
            if np.max(np.abs(ct.getDependentColumn(lab).to_numpy())) > 0.9:
                n_near += 1
    print(f"residuals near their +-1 bound (>0.9): {n_near}")

    mean_a = float(np.mean(wa))
    print("\nverdict on activation chatter:")
    if mean_a < ALT_WARN:
        print("  OK - activations are smooth; muscle forces are trustworthy.")
    elif mean_a < ALT_BAD:
        print("  WIGGLY - some mesh-scale ripple. Refine the mesh "
              "(~2 ms) before drawing muscle-level conclusions.")
    else:
        print("  CHATTERING - the mesh is too coarse for the activation "
              "dynamics (need h << 10 ms). Do NOT use these muscle "
              "activations/forces; re-run at ~2 ms (or solve in sub-windows).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
