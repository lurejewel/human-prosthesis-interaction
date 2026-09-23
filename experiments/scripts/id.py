"""Inverse Dynamics on the SQR level-walking trial -> joint moments.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\id.py"
    ... --start 170 --end 240                  # defaults

    level_walking_ik_filtered.mot  ->  ResultsInverseDynamics/level_walking_id.sto

Consumes exactly two band-limited inputs, and they must stay band-limited the
same way or the residual does not cancel:

    kinematics      level_walking_ik_filtered.mot   (ik_filter.py, 15 Hz)
    external loads  level_walking_grf_filtered.xml  (grf_filter.py, 15 Hz)

The tool applies a third 15 Hz filter of its own, configured in
Setup_InverseDynamics.xml (forces_to_exclude = Muscles).  A raw IK file must
never be passed here: inverse dynamics differentiates the coordinates twice, so
unfiltered noise becomes enormous joint accelerations.

Behaviour notes:
    - A previous result file is copied to <file>.bak before it is overwritten.
    - The <results_directory> / output name come from Setup_InverseDynamics.xml;
      `output` below is used for the existence check, the backup and the sanity
      report, and must match what the setup writes.
    - The sanity report prints peak pelvis residual force and peak joint moments
      inside the window, as a quick "does this look like walking" check.

Library use:
    from id import run_id
    res = run_id(coordinates=my_filtered_ik, t0=170.0, t1=240.0)
"""
import argparse
import os
import sys
from collections import namedtuple

import numpy as np

import opensim as osim

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import backup_if_exists, read_table, column_names   # noqa: E402
from ik_filter import IK_FILT_MOT                                   # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

MODEL = os.path.join(DATA_DIR, "SQR_simbody.osim")
ID_SETUP = os.path.join(DATA_DIR, "Setup_InverseDynamics.xml")
GRF_FILT_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")

ID_RESULTS = os.path.join(DATA_DIR, "ResultsInverseDynamics")
ID_OUT = os.path.join(ID_RESULTS, "level_walking_id.sto")

DEFAULT_TIME_RANGE = (170.0, 240.0)     # s, analysis window

# Reported by the sanity check: pelvis residual and the six sagittal moments.
SANITY_FORCES = ["pelvis_tx_force", "pelvis_ty_force"]
SANITY_MOMENTS = ["hip_flexion_r_moment", "hip_flexion_l_moment",
                  "knee_extension_r_moment", "knee_extension_l_moment",
                  "ankle_dorsiflexion_r_moment", "ankle_dorsiflexion_l_moment"]

IdResult = namedtuple("IdResult", ["path", "rows", "columns"])


# ---------------------------------------------------------------------------
# The tool
# ---------------------------------------------------------------------------
def run_id(coordinates=IK_FILT_MOT, output=ID_OUT, t0=DEFAULT_TIME_RANGE[0],
           t1=DEFAULT_TIME_RANGE[1], model=MODEL, external_loads=GRF_FILT_XML,
           setup=ID_SETUP, results_dir=ID_RESULTS, backup=True, sanity=True,
           verbose=True):
    """Run InverseDynamicsTool over [t0, t1]; write the .sto into `results_dir`.

    `coordinates` must be the FILTERED kinematics (see the module docstring).
    Returns an IdResult.
    """
    coordinates = os.path.abspath(coordinates)
    output = os.path.abspath(output)
    results_dir = os.path.abspath(results_dir)

    if not os.path.exists(coordinates):
        sys.exit("filtered kinematics not found: {}\n"
                 "Run ik.py + ik_filter.py (or run_ik_id.py) first.".format(coordinates))
    for label, path in (("model", model), ("setup", setup),
                        ("external loads", external_loads)):
        if not os.path.exists(path):
            sys.exit("ID {} file not found: {}".format(label, path))
    if backup:
        bak = backup_if_exists(output)
        if bak and verbose:
            print("   backed up existing result -> {}".format(os.path.basename(bak)))

    os.makedirs(results_dir, exist_ok=True)
    tool = osim.InverseDynamicsTool(setup)
    tool.setModelFileName(model)
    tool.setCoordinatesFileName(coordinates)
    tool.setExternalLoadsFileName(external_loads)
    tool.setResultsDir(results_dir)
    tool.setStartTime(t0)
    tool.setEndTime(t1)
    tool.run()

    if not os.path.exists(output):
        raise RuntimeError("ID did not produce output: " + output)

    _, col_labels, data = read_table(output)
    if verbose:
        print("   ID output: {} rows x {} cols".format(data.shape[0], data.shape[1]))
    if sanity:
        report_sanity(output, t0, t1)
    return IdResult(output, int(data.shape[0]), data.shape[1])


def report_sanity(path, t0, t1, verbose=True):
    """Peak pelvis residual force and peak joint moments inside the window."""
    _, col_labels, data = read_table(path)
    cols = column_names(col_labels)
    if len(cols) != data.shape[1]:
        cols = ["time"] + ["col{}".format(i + 1) for i in range(data.shape[1] - 1)]
    idx = {name: i for i, name in enumerate(cols)}
    mask = (data[:, 0] >= t0) & (data[:, 0] <= t1)

    def maxabs_in_window(name):
        return float(np.max(np.abs(data[mask, idx[name]])))

    if not verbose:
        return
    print("   [sanity] window [{:.0f} {:.0f}] s:".format(t0, t1))
    for j in SANITY_FORCES:
        if j in idx:
            print("     max|{}| = {:.1f} N".format(j, maxabs_in_window(j)))
    for j in SANITY_MOMENTS:
        if j in idx:
            print("     max|{}| = {:.1f} Nm".format(j, maxabs_in_window(j)))


# ---------------------------------------------------------------------------
# Main -- run this stage on its own
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Inverse Dynamics (filtered kinematics + filtered GRF)")
    ap.add_argument("--coordinates", default=IK_FILT_MOT,
                    help="FILTERED IK .mot to read (default: level_walking_ik_filtered.mot)")
    ap.add_argument("--output", default=ID_OUT, help="result .sto to write")
    ap.add_argument("--start", type=float, default=DEFAULT_TIME_RANGE[0],
                    help="start time in s (default {})".format(DEFAULT_TIME_RANGE[0]))
    ap.add_argument("--end", type=float, default=DEFAULT_TIME_RANGE[1],
                    help="end time in s (default {})".format(DEFAULT_TIME_RANGE[1]))
    args = ap.parse_args()

    print("Inverse Dynamics")
    print("  model       :", MODEL)
    print("  kinematics  :", args.coordinates)
    print("  external    :", GRF_FILT_XML)
    print("  setup       :", ID_SETUP)
    print("  window      : {:.0f}-{:.0f} s".format(args.start, args.end))

    res = run_id(args.coordinates, args.output, args.start, args.end)
    print("\nWrote {}.".format(res.path))


if __name__ == "__main__":
    main()
