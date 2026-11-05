"""Remove the heel-strike shear-force artifact from the measured GRF.

At heel strike the force plate rings: the anterior-posterior shear (vx) and the
medio-lateral shear (vz) briefly point in the *wrong* direction, then reverse.
This is a measurement artifact, not a real friction force, and it produces
large, simultaneous muscle-excitation spikes in MocoInverse.

Rule used here (early stance = the first `--early-window` seconds after contact):

    right foot:  vx <= 0  and  vz <= 0
    left  foot:  vx <= 0  and  vz >= 0

Wrong-sign runs are not hard-clipped to zero (that leaves a kink).  Each run is
replaced by a monotone PCHIP transition between the last valid sample before
the run and the first valid sample after it, so the corrected signal is smooth
and keeps the expected sign.  Everything outside the early-stance windows is
left untouched.

Outputs
-------
level_walking_grf_fixed.mot   corrected forces (same header/labels)
level_walking_grf_fixed.xml   ExternalLoads XML pointing at the fixed .mot
grf_heel_strike_fix_report.csv  every corrected interval

Run with the `opensim_scripting` conda environment (numpy/scipy are enough;
OpenSim itself is not needed):

    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe \
        experiments\\scripts\\fix_grf_heel_strike.py
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import read_table, write_table, write_loads_xml  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

DEFAULT_IN_MOT = os.path.join(DATA_DIR, "level_walking_grf_filtered.mot")
DEFAULT_IN_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")
DEFAULT_OUT_MOT = os.path.join(DATA_DIR, "level_walking_grf_fixed.mot")
DEFAULT_OUT_XML = os.path.join(DATA_DIR, "level_walking_grf_fixed.xml")


def contiguous_runs(mask):
    """[(first, last), ...] for each contiguous True run in a bool array."""
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        return []
    parts = np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1)
    return [(int(p[0]), int(p[-1])) for p in parts]


def heel_strikes(vy, t, threshold, min_gap):
    """Indices where vy rises above `threshold`, debounced by `min_gap`."""
    above = vy > threshold
    rising = np.flatnonzero((~above[:-1]) & above[1:]) + 1
    hs = []
    for k in rising:
        k = int(k)
        if not hs or t[k] - t[hs[-1]] > min_gap:
            hs.append(k)
    if above[0] and (not hs or t[0] < t[hs[0]] - min_gap):
        hs.insert(0, 0)
    return hs


def fix_component(x, t, valid, i0, i1, fs, rows, side, comp):
    """Replace wrong-sign runs in x[i0:i1] by smooth sign-preserving PCHIP."""
    x_out = x.copy()
    # include the sample just before the threshold crossing, because the
    # artifact often starts while vy is still below the detection threshold
    j0 = max(0, i0 - 1)
    bad = np.zeros(x.size, dtype=bool)
    bad[j0:i1 + 1] = ~valid(x[j0:i1 + 1])

    for a, b in contiguous_runs(bad):
        orig = x[a:b + 1].copy()

        # left anchor: the last sample before the run; 0 if it is not valid
        li = a - 1
        if li < 0:
            li = 0
            lv = 0.0
        else:
            lv = float(x[li])
            if not bool(valid(np.array([lv]))[0]):
                lv = 0.0

        # right anchor: the first valid sample after the run (search up to
        # 0.1 s past the early-stance window, in case the run ends at i1)
        search_end = min(x.size - 1, max(i1, b) + int(round(0.1 * fs)))
        ri = b + 1
        rv = None
        while ri <= search_end:
            if bool(valid(np.array([x[ri]]))[0]):
                rv = float(x[ri])
                break
            ri += 1
        if rv is None:
            ri = min(x.size - 1, b + 1)
            rv = 0.0

        if ri <= li or t[ri] <= t[li]:
            x_out[a:b + 1] = 0.0
        else:
            p = PchipInterpolator(np.array([t[li], t[ri]]),
                                  np.array([lv, rv]))
            vals = p(t[a:b + 1])
            # the anchors have the right sign, so PCHIP should too; guard anyway
            vals = np.where(valid(vals), vals, 0.0)
            x_out[a:b + 1] = vals

        rows.append(dict(side=side, component=comp,
                         hs_time=float(t[i0]),
                         t_start=float(t[a]), t_end=float(t[b]),
                         n_samples=int(b - a + 1),
                         orig_min=float(orig.min()),
                         orig_max=float(orig.max())))

    return x_out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fix the heel-strike vx/vz sign artifact in the GRF")
    ap.add_argument("--input", default=DEFAULT_IN_MOT)
    ap.add_argument("--xml", default=DEFAULT_IN_XML)
    ap.add_argument("--output", default=DEFAULT_OUT_MOT)
    ap.add_argument("--out-xml", default=DEFAULT_OUT_XML)
    ap.add_argument("--early-window", type=float, default=0.2,
                    help="seconds after heel strike to enforce the sign rule")
    ap.add_argument("--hs-threshold", type=float, default=2.0,
                    help="vertical GRF [N] used to detect initial contact")
    ap.add_argument("--min-gap", type=float, default=0.3,
                    help="minimum time between detected heel strikes [s]")
    args = ap.parse_args(argv)

    header, labels, data = read_table(args.input)
    t = data[:, 0].copy()
    fs = 1.0 / float(np.median(np.diff(t)))
    lab = [c.strip() for c in labels.split()]
    print(f"input  : {args.input}")
    print(f"rows   : {data.shape[0]}  fs = {fs:.1f} Hz")
    print(f"rule   : early {args.early_window:g} s; "
          f"right vx<=0,vz<=0; left vx<=0,vz>=0")

    rows = []
    for side, prefix, vz_valid in [
            ("left", "ground_force1", lambda v: v >= 0),
            ("right", "ground_force2", lambda v: v <= 0)]:
        ix = lab.index(prefix + "_vx")
        iy = lab.index(prefix + "_vy")
        iz = lab.index(prefix + "_vz")
        vy = data[:, iy]
        hs = heel_strikes(vy, t, args.hs_threshold, args.min_gap)
        print(f"\n{side}: {len(hs)} heel strikes "
              f"({', '.join(f'{t[k]:.3f}' for k in hs)})")

        for k0 in hs:
            i1 = int(np.searchsorted(t, t[k0] + args.early_window,
                                     side="right") - 1)
            i1 = min(i1, len(t) - 1)
            data[:, ix] = fix_component(data[:, ix], t,
                                        lambda v: v <= 0, k0, i1, fs,
                                        rows, side, "vx")
            data[:, iz] = fix_component(data[:, iz], t, vz_valid,
                                        k0, i1, fs, rows, side, "vz")

    write_table(args.output, header, labels, data)
    write_loads_xml(args.out_xml, args.output, template_xml=args.xml)
    print(f"\nwrote  : {args.output}")
    print(f"wrote  : {args.out_xml}")

    rep = os.path.join(os.path.dirname(os.path.abspath(args.output)),
                       "grf_heel_strike_fix_report.csv")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("side,component,hs_time,t_start,t_end,n_samples,orig_min,orig_max\n")
        for r in rows:
            f.write("{side},{component},{hs_time:.6f},{t_start:.6f},"
                    "{t_end:.6f},{n_samples},{orig_min:.4f},{orig_max:.4f}\n"
                    .format(**r))
    print(f"report : {rep}")
    print(f"corrected intervals: {len(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
