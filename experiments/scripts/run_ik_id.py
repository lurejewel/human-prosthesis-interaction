"""IK + ID on the SQR level-walking trial -- pipeline wiring only.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\run_ik_id.py"

Prerequisite, the filtered external loads:
    python "experiments\\scripts\\grf_filter.py"
which writes level_walking_grf_filtered.mot / .xml.  This pipeline does NOT
filter the GRF; it consumes that file so that it and run_ik_rra.py always see
identical external loads.

This file only names the inputs/outputs and calls the stages; all the work lives
in the modules next to it:

    [1/3]  ik.run_ik()        IK          -> level_walking_ik.mot          (RAW)
    [2/3]  ik_filter.filter_ik()  15 Hz   -> level_walking_ik_filtered.mot
    [3/3]  id.run_id()        ID          -> ResultsInverseDynamics/level_walking_id.sto
                                              (reads the FILTERED kinematics)

There are no options: every run redoes all three stages, so the outputs cannot
mix old and new stages.  Each stage copies the file it overwrites to .bak first.

The 15 Hz cutoff has to agree in three places, otherwise the external loads and
the twice-differentiated kinematics are band-limited differently and the
inverse-dynamics residual does not cancel:
    grf_filter.GRF_FILTER_FC            (the GRF itself)
    ik_filter.KIN_FILTER_FC             (the kinematics, stage 2)
    Setup_InverseDynamics.xml           (the tool's own coordinate filter)

Each stage is also runnable on its own: ik.py, ik_filter.py, id.py.
"""
import os
import sys

import opensim as osim

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import GRF_FILT_MOT, GRF_FILT_XML                 # noqa: E402
from ik import DATA_DIR, IK_RAW_MOT, run_ik                       # noqa: E402
from ik_filter import IK_FILT_MOT, KIN_FILTER_FC, filter_ik       # noqa: E402
from id import ID_OUT, run_id                                     # noqa: E402

# This pipeline's own inputs and outputs.
IK_RAW = IK_RAW_MOT                     # raw IK       (written by stage 1)
IK_FILT = IK_FILT_MOT                   # filtered IK  (written by stage 2)

TIME_RANGE = (170.0, 240.0)             # s, analysis window


def check_inputs():
    """The only input this pipeline needs from outside: the filtered GRF."""
    missing = [p for p in (GRF_FILT_MOT, GRF_FILT_XML) if not os.path.exists(p)]
    if missing:
        sys.exit("Missing input files:\n  " + "\n  ".join(missing) +
                 "\n\nRun grf_filter.py first.")


def main():
    print("OpenSim version :", osim.GetVersion())
    print("Script          :", os.path.abspath(__file__))
    print("Data dir        :", DATA_DIR)
    print("Window          : {:.0f}-{:.0f} s".format(*TIME_RANGE))
    print("Cutoff          : {:.1f} Hz (kinematics; GRF comes from grf_filter.py)"
          .format(KIN_FILTER_FC))
    check_inputs()

    # OpenSim resolves any relative path inside the setup XMLs against the cwd,
    # so work from the data directory exactly as the GUI does
    os.chdir(DATA_DIR)

    print("\n[1/3] Inverse Kinematics ...")
    run_ik(IK_RAW, TIME_RANGE[0], TIME_RANGE[1])
    print("   raw kinematics -> {}".format(IK_RAW))

    print("\n[2/3] Low-pass filtering IK coordinates @ {:.1f} Hz ...".format(KIN_FILTER_FC))
    filter_ik(IK_RAW, IK_FILT, cutoff=KIN_FILTER_FC, backup=True)
    print("   filtered kinematics -> {}".format(IK_FILT))

    print("\n[3/3] Inverse Dynamics (filtered kinematics + filtered GRF) ...")
    run_id(IK_FILT, ID_OUT, TIME_RANGE[0], TIME_RANGE[1])

    print("\nPipeline finished.")
    print("  raw IK      : {}".format(IK_RAW))
    print("  filtered IK : {}  <- use this one".format(IK_FILT))
    print("  ID result   : {}".format(ID_OUT))


if __name__ == "__main__":
    main()
