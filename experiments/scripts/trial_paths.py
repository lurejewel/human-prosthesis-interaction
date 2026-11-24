"""Resolve the trial files that live next to this script.

Copied verbatim into every trial folder; everything is derived from the
folder the script sits in, so no per-group editing is needed.
"""

import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir, os.pardir))
TRIAL = os.path.basename(HERE)


def one(pattern):
    hits = sorted(glob.glob(os.path.join(HERE, pattern)))
    if not hits:
        raise SystemExit(f"{TRIAL}: no file matching {pattern}")
    return hits[0]


GRF_RAW_MOT = one("*_grf.mot")
GRF_RAW_XML = one("*_grf.xml")
IK_RAW_MOT = one("*_ik.mot")
MODEL_FILE = one("*_simbody.osim")
STEM = os.path.basename(GRF_RAW_MOT)[:-len("_grf.mot")]

GRF_FILT_MOT = os.path.join(HERE, STEM + "_grf_filtered.mot")
GRF_FILT_XML = os.path.join(HERE, STEM + "_grf_filtered.xml")
GRF_FIXED_MOT = os.path.join(HERE, STEM + "_grf_fixed.mot")
GRF_FIXED_XML = os.path.join(HERE, STEM + "_grf_fixed.xml")
IK_FILT_MOT = os.path.join(HERE, STEM + "_ik_filtered.mot")
SEMG_PROCESSED = os.path.join(HERE, STEM + "_semg_processed.csv")
EMG_REF_STO = os.path.join(HERE, STEM + "_emg_reference.sto")
EMG_REF_CSV = os.path.join(HERE, STEM + "_emg_reference.csv")
OUT_DIR = os.path.join(HERE, "moco_outputs")


def window():
    path = os.path.join(HERE, "steady_state_summary.json")
    if not os.path.isfile(path):
        raise SystemExit(f"{TRIAL}: run detect_steady_state.py first")
    return json.load(open(path))["analysis_window"]