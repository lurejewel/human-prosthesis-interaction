"""IK + RRA on the SQR level-walking trial -- pipeline wiring only.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\scripts\\run_ik_rra.py"
    ... --passes 5 --start 170 --end 240          # defaults

Prerequisite (this pipeline does NOT filter the GRF):
    python "experiments\\scripts\\grf_filter.py"      # -> level_walking_grf_filtered.mot/.xml

This file only names the inputs/outputs and calls the stages; all the work lives
in the modules next to it:

    [1/4]  ik.run_ik()                   IK      -> level_walking_ik.mot          (RAW)
    [2/4]  ik_filter.filter_ik()         15 Hz   -> level_walking_ik_filtered.mot
    [3/4]  rra.build_model()                      -> RRA_setup/ (model + setup XMLs)
           rra.write_window_files()               -> RRA_setup/ (window files)
    [4/4]  rra.run_rra_passes()          RRA     -> ResultsRRA/RRA_pass<k>_*

Stages 1-2 are the very same IK and the same 15 Hz filter as run_ik_id.py -- the
modules are shared -- so this pipeline writes those two files rather than a
second, byte-identical copy under a _rra name.  Both are overwritten on every
run; with a non-default --start/--end the shared files then hold THIS window's
kinematics, so re-run run_ik_id.py afterwards if you need its 170-240 s versions
back.

Options mirror the RRA stage (see rra.py for what each parameter does):
    --passes K    number of RRA rounds (default 5); 0 stops after preparing
    --start/--end solve window in s (default 170-240)
    --resume K    start at pass K from ResultsRRA/RRA_pass<K-1>_adjusted.osim

RRA is fed the FILTERED kinematics and applies its own 15 Hz lowpass on top
(setLowpassCutoffFrequency), matching the GRF and Setup_InverseDynamics.xml.
"""
import argparse
import os
import sys

import opensim as osim

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ik                                                             # noqa: E402
import rra                                                            # noqa: E402
from grf_filter import GRF_FILT_MOT, GRF_FILT_XML                     # noqa: E402
from ik import IK_RAW_MOT                                             # noqa: E402
from ik_filter import IK_FILT_MOT, KIN_FILTER_FC, filter_ik           # noqa: E402

# The shared kinematics, written by stages 1-2 (see the docstring).
IK_RAW = IK_RAW_MOT                     # level_walking_ik.mot, raw
IK_FILT = IK_FILT_MOT                   # level_walking_ik_filtered.mot, 15 Hz


def check_inputs():
    """The only input this pipeline needs from outside: the filtered GRF."""
    required = [GRF_FILT_MOT, GRF_FILT_XML, ik.MODEL, ik.MARKERS, ik.IK_SETUP]
    missing = [p for p in required if not os.path.exists(p)]
    if missing:
        sys.exit("Missing input files:\n  " + "\n  ".join(missing) +
                 "\n\nRun grf_filter.py first if the filtered GRF is missing.")


def main():
    ap = argparse.ArgumentParser(description="IK + RRA on the SQR walking trial")
    ap.add_argument("--passes", type=int, default=rra.DEFAULT_PASSES,
                    help="number of RRA passes (default {})".format(rra.DEFAULT_PASSES))
    ap.add_argument("--start", type=float, default=rra.DEFAULT_TIME_RANGE[0])
    ap.add_argument("--end", type=float, default=rra.DEFAULT_TIME_RANGE[1])
    ap.add_argument("--resume", type=int, default=1, metavar="K",
                    help="start at pass K, reusing ResultsRRA/RRA_pass<K-1>_adjusted.osim")
    args = ap.parse_args()
    t0, t1 = args.start, args.end

    print("OpenSim version :", osim.GetVersion())
    print("Script          :", os.path.abspath(__file__))
    print("Window          : {:.0f}-{:.0f} s".format(t0, t1))
    print("RRA passes      :", args.passes)
    print("Cutoff          : {:.1f} Hz".format(KIN_FILTER_FC))
    check_inputs()

    print("\n[1/4] Inverse Kinematics ...")
    ik.run_ik(IK_RAW, t0, t1)
    print("   raw kinematics -> {}".format(IK_RAW))

    print("\n[2/4] Low-pass filtering IK coordinates @ {:.1f} Hz ...".format(KIN_FILTER_FC))
    filter_ik(IK_RAW, IK_FILT, cutoff=KIN_FILTER_FC, backup=True)
    print("   filtered kinematics -> {}".format(IK_FILT))

    print("\n[3/4] Building the RRA model, setup files and window files ...")
    rra.build_model()
    _, grf_window_xml, ik_window = rra.write_window_files(IK_FILT, t0, t1)

    print("\n[4/4] Running RRA x{} ...".format(args.passes))
    history = rra.run_rra_passes(t0, t1, n_passes=args.passes,
                                 ik_window=ik_window,
                                 grf_window_xml=grf_window_xml,
                                 resume=args.resume)
    rra.print_summary(history, args.passes)


if __name__ == "__main__":
    main()
