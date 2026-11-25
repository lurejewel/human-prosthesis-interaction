"""Run the MocoInverse pre-processing chain for the trial folder this file sits in.

Stages (all inputs and outputs stay inside this folder):
  grf   raw GRF -> 15 Hz zero-phase, stance-segmented, debounced, tapered
  ik    raw IK -> 15 Hz zero-phase (same cutoff as the GRF)
  semg  raw sEMG -> despiked / 10 Hz HP / rectified / 5 Hz LP / max-normalised
  emg   envelope -> 100 Hz MocoControlTrackingGoal reference (.sto)

B-form inner/outer excitation estimation is run separately with
`run_bform_outer.py`; this script never calls MocoInverse directly.

`fix` (heel-strike shear-sign correction) is deliberately NOT in the default
chain: its sign rule encodes SQR's walking direction (+x).  The JW trials walk
in -x, so their early-stance shear already has the correct sign and applying
the rule would invert ~600 valid intervals.  Run it only via --stages fix if a
genuine wrong-sign artifact is ever identified.

Usage:
    python run_pipeline.py                 # pre-processing only
    python run_bform_outer.py --workers 4  # B-form inner/outer excitation solve
"""

import argparse
import os
import sys

import trial_paths as T

sys.path.insert(0, T.HERE)
os.chdir(T.HERE)          # opensim.log / solver sentinel stay in this folder

FC = 15.0                 # shared GRF / IK cutoff, Hz
ORDER = 4


def stage_grf():
    import grf_filter as gf
    gf.GRF_RAW_MOT, gf.GRF_RAW_XML = T.GRF_RAW_MOT, T.GRF_RAW_XML
    gf.GRF_FILT_MOT, gf.GRF_FILT_XML = T.GRF_FILT_MOT, T.GRF_FILT_XML
    gf.TIME_RANGE = tuple(T.window())
    # template_xml is a def-time default, so overriding the module global is not enough
    gf.write_loads_xml.__defaults__ = (T.GRF_RAW_XML,)
    gf.filter_grf(quiet=True)


def stage_fix():
    import fix_grf_heel_strike as fx
    fx.main(["--input", T.GRF_FILT_MOT, "--xml", T.GRF_FILT_XML,
             "--output", T.GRF_FIXED_MOT, "--out-xml", T.GRF_FIXED_XML])


def stage_ik():
    import ik_filter as ikf
    ikf.filter_ik(src=T.IK_RAW_MOT, dst=T.IK_FILT_MOT,
                  cutoff=FC, order=ORDER, quiet=True)


def stage_semg():
    import process_semg_jw as ps
    ps.main(["--data-dir", T.HERE])


def stage_emg():
    import make_emg_reference as mer
    a0, a1 = T.window()
    mer.main(["--input", T.SEMG_PROCESSED,
              "--start", f"{a0 - 1.0:.3f}", "--end", f"{a1 + 1.0:.3f}",
              "--dt", "0.01", "--out", T.EMG_REF_STO, "--out-csv", T.EMG_REF_CSV])


STAGES = [("grf", stage_grf), ("fix", stage_fix), ("ik", stage_ik),
          ("semg", stage_semg), ("emg", stage_emg)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stages", default="grf,ik,semg,emg",
                    help="comma list of pre-processing stages to run")
    args = ap.parse_args()
    wanted = [s.strip() for s in args.stages.split(",") if s.strip()]
    for name, fn in STAGES:
        if name in wanted:
            print(f"\n=== [{T.TRIAL}] stage {name} ===", flush=True)
            fn()
    print(f"\n[{T.TRIAL}] done. window={T.window()}")


if __name__ == "__main__":
    main()