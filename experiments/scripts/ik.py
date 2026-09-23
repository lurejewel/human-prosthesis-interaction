"""Inverse Kinematics on the SQR level-walking trial -> RAW coordinate .mot.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\ik.py"
    ... --start 170 --end 240                  # defaults

Writes level_walking_ik.mot: the RAW InverseKinematicsTool result.  Nothing may
differentiate it directly -- send it through ik_filter.py first (both pipelines
do that automatically).  See ik_filter.py on the naming:

    <trial>_ik.mot            raw IK       <- this module
    <trial>_ik_filtered.mot   15 Hz filtered, what ID / RRA / Moco read

Both pipelines call run_ik() with the same model, markers and setup, so the tool
configuration lives here exactly once and both write the same file:
    run_ik_id.py    -> level_walking_ik.mot
    run_ik_rra.py   -> level_walking_ik.mot      (same file, no _rra copy)
Only the time window may differ (run_ik_rra.py --start/--end).

Inputs (all in experiments/data/SQR_walking/):
    SQR_simbody.osim    model
    level_walking.trc   measured marker trajectories
    Setup_IK.xml        marker weights, coordinate limits, solver settings

Behaviour notes:
    - Every file is handed to the tool as an absolute path, so the caller does
      not have to chdir into the data directory first.
    - A previous output file is copied to <file>.bak before it is overwritten.
    - OpenSim writes the per-marker error report into the *current working
      directory*, not next to the model, so it is looked for in both the output
      directory and the cwd; the newest match wins.
    - Fit quality on this trial: marker error RMS ~0.03 m, max ~0.06 m.

Library use:
    from ik import run_ik
    res = run_ik(output=".../level_walking_ik.mot", t0=170.0, t1=240.0)
    print(res.rows, res.marker_errors)
"""
import argparse
import glob
import os
import sys
from collections import namedtuple

import opensim as osim

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import backup_if_exists, read_table               # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

MODEL = os.path.join(DATA_DIR, "SQR_simbody.osim")
MARKERS = os.path.join(DATA_DIR, "level_walking.trc")
IK_SETUP = os.path.join(DATA_DIR, "Setup_IK.xml")

IK_RAW_MOT = os.path.join(DATA_DIR, "level_walking_ik.mot")   # raw IK, do NOT differentiate
DEFAULT_TIME_RANGE = (170.0, 240.0)                           # s, analysis window

# What run_ik() hands back to the caller.
IkResult = namedtuple("IkResult", ["path", "rows", "t_start", "t_end",
                                   "marker_errors"])


# ---------------------------------------------------------------------------
# The tool
# ---------------------------------------------------------------------------
def _newest_marker_errors(paths):
    """Of `paths`, the most recently modified one, or None."""
    hits = [p for p in paths if os.path.exists(p)]
    if not hits:
        return None
    return max(hits, key=os.path.getmtime)


def find_marker_errors(output):
    """Locate the per-marker error report the tool just wrote.

    OpenSim drops it in the process cwd, which is not necessarily the output
    directory, so both are searched.
    """
    dirs = [os.path.dirname(os.path.abspath(output)), os.getcwd()]
    seen, cands = set(), []
    for d in dirs:
        d = os.path.abspath(d)
        if d in seen:
            continue
        seen.add(d)
        cands += glob.glob(os.path.join(d, "*marker_errors*.sto"))
    return _newest_marker_errors(cands)


def run_ik(output=IK_RAW_MOT, t0=DEFAULT_TIME_RANGE[0], t1=DEFAULT_TIME_RANGE[1],
           model=MODEL, markers=MARKERS, setup=IK_SETUP, backup=True, verbose=True):
    """Run InverseKinematicsTool over [t0, t1] and write the RAW result to `output`.

    Returns an IkResult.  Nothing is filtered here; see the module docstring.
    """
    output = os.path.abspath(output)
    for label, path in (("model", model), ("marker", markers), ("setup", setup)):
        if not os.path.exists(path):
            sys.exit("IK {} file not found: {}".format(label, path))
    if backup:
        bak = backup_if_exists(output)
        if bak and verbose:
            print("   backed up existing IK -> {}".format(os.path.basename(bak)))

    tool = osim.InverseKinematicsTool(setup)
    tool.set_model_file(model)
    tool.set_marker_file(markers)
    tool.set_output_motion_file(output)
    tool.setStartTime(t0)
    tool.setEndTime(t1)
    tool.run()

    if not os.path.exists(output):
        sys.exit("IK did not produce " + output)
    _, _, data = read_table(output)

    errors = find_marker_errors(output)
    if verbose:
        print("   IK output: {} rows, time {:.2f} - {:.2f} s  (raw)".format(
            data.shape[0], data[0, 0], data[-1, 0]))
        print("   marker error report: {}".format(errors or "(none found)"))
    return IkResult(output, int(data.shape[0]), float(data[0, 0]), float(data[-1, 0]),
                    errors)


# ---------------------------------------------------------------------------
# Main -- run this stage on its own
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Inverse Kinematics (raw output)")
    ap.add_argument("--output", default=IK_RAW_MOT,
                    help="raw IK .mot to write (default: the SQR level-walking IK)")
    ap.add_argument("--start", type=float, default=DEFAULT_TIME_RANGE[0],
                    help="start time in s (default {})".format(DEFAULT_TIME_RANGE[0]))
    ap.add_argument("--end", type=float, default=DEFAULT_TIME_RANGE[1],
                    help="end time in s (default {})".format(DEFAULT_TIME_RANGE[1]))
    args = ap.parse_args()

    print("Inverse Kinematics (raw)")
    print("  model   :", MODEL)
    print("  markers :", MARKERS)
    print("  setup   :", IK_SETUP)
    print("  window  : {:.0f}-{:.0f} s".format(args.start, args.end))

    res = run_ik(args.output, args.start, args.end)
    print("\nWrote {} (raw).".format(res.path))
    print("Filter it before use:  ik_filter.py --input \"{}\"".format(res.path))


if __name__ == "__main__":
    main()
