"""Unit test of the segmentation / stitching logic (no solver involved)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location(
    "moco_inverse", os.path.join(HERE, "moco_inverse.py"))
mi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mi)

# ---- 1. plan covers [t0, t1] exactly, with burn-in on the solve windows ----
plan = mi.build_segment_plan(180.0, 200.0, 2.0, 0.5, 180.0, 240.0)
print(f"segments: {len(plan)}")
print(f"  {'k':>3}{'retain':>20}{'solve':>24}{'len':>7}")
for p in plan:
    print(f"  {p['k']:>3}   [{p['retain'][0]:7.3f},{p['retain'][1]:7.3f}]"
          f"   [{p['solve'][0]:7.3f},{p['solve'][1]:7.3f}]"
          f"{p['solve'][1]-p['solve'][0]:>7.2f}")

assert all(p["retain"][0] == (180.0 + 2.0 * p["k"]) for p in plan)
assert all(p["retain"][1] - p["retain"][0] == 2.0 for p in plan)
assert plan[0]["retain"][0] == 180.0 and plan[-1]["retain"][1] == 200.0
for a, b in zip(plan[:-1], plan[1:]):
    assert a["retain"][1] == b["retain"][0], "gap or overlap in retained slices"
print("  retained slices tile [180,200] exactly: OK")

# burn-in really exists on both sides (except where the data clamps it)
assert plan[0]["solve"][0] == 180.0 and plan[0]["solve"][1] == 182.5
assert plan[1]["solve"][0] == 181.5 and plan[1]["solve"][1] == 184.5
assert plan[-1]["solve"][1] == 200.5
print("  burn-in windows: OK  (seg0 has no data before 180, seg9 none after 200)")

# degenerate / edge cases
assert len(mi.build_segment_plan(180, 180.3, 2.0, 0.5, 180, 240)) == 1
assert mi.build_segment_plan(180, 180.3, 2.0, 0.5, 180, 240)[0]["solve"] == (180, 180.3)
assert len(mi.build_segment_plan(180, 190, 0, 0.5, 180, 240)) == 1   # --segment 0
p = mi.build_segment_plan(180, 185, 2.0, 0.5, 180, 240)
assert len(p) == 3 and p[-1]["retain"] == (184.0, 185.0)             # last slice short
print("  edge cases (window <= segment, segment=0, non-multiple end): OK")

# ---- 2. ownership: every grid point belongs to exactly one segment --------
grid = np.arange(180.0, 200.0 + 0.0005, 0.001)
owner = mi.owners_for(grid, plan)
assert len(owner) == len(grid)
counts = np.bincount(owner, minlength=len(plan))
assert counts.sum() == len(grid)
assert owner[0] == 0 and owner[-1] == len(plan) - 1
print(f"\nownership over {len(grid)} grid points ({counts.tolist()}), "
      f"monotone={bool(np.all(np.diff(owner) >= 0))}")

# ---- 3. stitch() with synthetic segments --------------------------------
# each "segment" returns a smooth function plus a deliberate burn-in error at
# its start, so the seam diagnostic must detect it.
def fake_segment(k, ts, bump):
    a = np.sin(0.5 * (ts - 180.0))[:, None]
    if bump:                       # burn-in error decaying with tau = 40 ms
        a = a + 0.3 * np.exp(-(ts - ts[0]) / 0.04)[:, None]
    ctrl = 0.5 * np.cos(0.5 * (ts - 180.0))[:, None]
    return (ts, ["/forceset/x/activation"], a,
            ["/forceset/x"], ctrl)


data = []
for p in plan:
    ts = np.arange(p["solve"][0], p["solve"][1] + 1e-9, 0.001)
    data.append(fake_segment(p["k"], ts, bump=(p["k"] > 0)))

st_lab, states, ct_lab, ctrls, seams, ai, ei = mi.stitch(data, plan, grid)
assert states.shape == (len(grid), 1) and not np.isnan(states).any()
print(f"\nstitch(): no NaNs, {len(seams)} seams")
print(f"  {'cut':>8}{'n':>6}{'act RMS':>11}{'act max':>11}")
for s in seams:
    print(f"  {s['cut']:>8.1f}{s['n']:>6}{s['act_rms']:>11.5f}{s['act_max']:>11.5f}")

# with a 0.5 s burn-in and tau=40 ms, the recovered burn-in error at a seam is
# 0.3*exp(-0.5/0.04) ~ 1e-6, so the seam diagnostics must be tiny
worst = max(s["act_max"] for s in seams)
assert worst < 1e-3, f"seam diagnostic did not wash out the burn-in: {worst}"
print(f"  worst seam mismatch = {worst:.2e}  -> burn-in washes the error out: OK")

# and the stitched signal must equal the clean reference (no bump)
ref = np.sin(0.5 * (grid - 180.0))
err = np.max(np.abs(states[:, 0] - ref))
assert err < 1e-3, f"stitched series deviates from reference: {err}"
print(f"  max deviation from reference = {err:.2e}: OK")

print("\nALL SEGMENTATION / STITCHING TESTS PASSED")
