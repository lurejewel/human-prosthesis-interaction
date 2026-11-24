"""Copy the trial pipeline scripts into each JW_* data folder.

Every copied script resolves its inputs and outputs from the folder it lives
in (see trial_paths.py), so a group folder becomes a self-contained trial:
run `python run_pipeline.py` inside it and nothing is written elsewhere.

Usage:
    python experiments/scripts/setup_group_pipeline.py            # all JW_* groups
    python experiments/scripts/setup_group_pipeline.py JW_run     # one group
"""

import glob
import os
import shutil
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MOCO_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "moco"))
DATA_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))

FROM_SCRIPTS = ["grf_filter.py", "fix_grf_heel_strike.py", "ik_filter.py",
                "make_emg_reference.py", "detect_steady_state.py",
                "process_semg_jw.py", "derotate_grf.py", "trial_paths.py",
                "run_pipeline.py", "visualize_preprocessing.m"]
FROM_MOCO = ["moco_inverse.py", "moco_bootstrap.py", "check_solution_quality.py"]


def main():
    targets = sys.argv[1:]
    if targets:
        folders = [os.path.join(DATA_DIR, t) for t in targets]
    else:
        folders = sorted(d for d in glob.glob(os.path.join(DATA_DIR, "JW_*"))
                         if os.path.isdir(d))
    for folder in folders:
        if not os.path.isdir(folder):
            sys.exit(f"not a directory: {folder}")
        for name in FROM_SCRIPTS:
            shutil.copy2(os.path.join(SCRIPT_DIR, name), os.path.join(folder, name))
        for name in FROM_MOCO:
            shutil.copy2(os.path.join(MOCO_DIR, name), os.path.join(folder, name))
        print(f"{os.path.basename(folder)}: copied {len(FROM_SCRIPTS) + len(FROM_MOCO)} scripts")


if __name__ == "__main__":
    main()