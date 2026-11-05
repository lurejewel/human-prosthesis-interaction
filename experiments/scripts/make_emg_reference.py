"""Build a MocoControlTrackingGoal reference table from the processed sEMG.

The processed sEMG file (`level_walking_semg_processed.csv`, produced by
`process_semg.m`) stores the rectified / 5 Hz low-passed / max-normalized
envelope in columns EMG7..EMG16.  This script resamples those 10 channels onto
a uniform grid and writes an OpenSim `TimeSeriesTable` whose column labels are
the **control names** used by MocoInverse:

    /forceset/rect_fem_r   (EMG7)
    /forceset/vasti_r      (EMG8)
    /forceset/hamstrings_r (EMG9)   <- biceps femoris sEMG, model "hamstrings"
    /forceset/tib_ant_r    (EMG10)
    /forceset/gastroc_r    (EMG11)
    ... and the same five for the left side (EMG12..EMG16)

Only these 10 muscles get a reference; the other 8 model muscles have no sEMG
and are left to the excitation-effort term.

Usage:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe \
        experiments\\scripts\\make_emg_reference.py \
        --start 199.5 --end 204.0 --dt 0.01
"""
from __future__ import annotations

import argparse
import os

import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")
SEMG_CSV = os.path.join(SCRIPT_DIR, "..", "moco", "outputs", "semg",
                        "level_walking_semg_processed.csv")

# (CSV column, model muscle / control name)
CHANNELS = [
    ("EMG7",  "rect_fem_r"),   # right rectus femoris
    ("EMG8",  "vasti_r"),      # right vastus lateralis
    ("EMG9",  "hamstrings_r"), # right biceps femoris -> model hamstrings
    ("EMG10", "tib_ant_r"),    # right tibialis anterior
    ("EMG11", "gastroc_r"),    # right gastrocnemius
    ("EMG12", "rect_fem_l"),
    ("EMG13", "vasti_l"),
    ("EMG14", "hamstrings_l"),
    ("EMG15", "tib_ant_l"),
    ("EMG16", "gastroc_l"),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=SEMG_CSV)
    ap.add_argument("--start", type=float, default=199.5,
                    help="reference start [s] (include the burn-in)")
    ap.add_argument("--end", type=float, default=204.0)
    ap.add_argument("--dt", type=float, default=0.01,
                    help="reference sampling interval [s] (100 Hz)")
    ap.add_argument("--out", default=os.path.join(
        SCRIPT_DIR, "..", "moco", "outputs",
        "emg_excitation_199p5_204p0.sto"))
    ap.add_argument("--out-csv", default=os.path.join(
        SCRIPT_DIR, "..", "moco", "outputs",
        "emg_excitation_199p5_204p0.csv"))
    args = ap.parse_args(argv)

    raw = np.loadtxt(args.input, delimiter=",", skiprows=1)
    t = raw[:, 0]
    header = open(args.input, "r", encoding="utf-8").readline().strip().split(",")
    col = {name: i for i, name in enumerate(header)}
    print(f"input : {args.input}")
    print(f"rows  : {len(t)}  span {t[0]:.3f}..{t[-1]:.3f} s")

    grid = np.arange(args.start, args.end + 1e-12, args.dt)
    data = np.empty((len(grid), len(CHANNELS) + 1))
    data[:, 0] = grid
    for j, (ch, muscle) in enumerate(CHANNELS, start=1):
        v = raw[:, col[ch]]
        data[:, j] = np.interp(grid, t, v)
        print(f"  {ch:5s} -> {muscle:14s}  "
              f"min={data[:, j].min():.3f} max={data[:, j].max():.3f}")

    labels = ["/forceset/" + m for _, m in CHANNELS]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write("version=1\n")
        f.write(f"nRows={len(grid)}\n")
        f.write(f"nColumns={len(labels) + 1}\n")
        f.write("inDegrees=no\nendheader\n")
        f.write("time\t" + "\t".join(labels) + "\n")
        for row in data:
            f.write("\t".join(f"{v:.8f}" for v in row) + "\n")
    print(f"wrote : {args.out}")

    with open(args.out_csv, "w", encoding="utf-8") as f:
        f.write("time," + ",".join(m for _, m in CHANNELS) + "\n")
        for row in data:
            f.write(",".join(f"{v:.8f}" for v in row) + "\n")
    print(f"wrote : {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
