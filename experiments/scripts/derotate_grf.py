"""DIAGNOSTIC ONLY: report the treadmill-frame rotation of a slope trial's GRF.

On a treadmill the centre of mass is stationary on average, so the mean total
GRF must point along true vertical.  The JW slope trials show a mean tilted
~15 deg about z with |mean| == BW.  This is EXPECTED and must NOT be corrected:
the world origin sits on the inclined treadmill, so markers and forces are both
expressed in the treadmill-aligned frame, and each group's JW_simbody.osim
already carries the matching rotated gravity (uphill [-2.5381, -9.4725, 0],
downhill [+2.5381, -9.4725, 0]).  Dynamics are therefore self-consistent as-is.

This script only writes grf_frame_report.json documenting the angle; it never
modifies force data.  Do not de-rotate the GRF.
"""

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trial_paths as T
from grf_filter import read_table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)

    header, labels, data = read_table(T.GRF_RAW_MOT)
    f = data[:, 1:4] + data[:, 10:13]
    mu = f.mean(0)
    theta = float(np.arctan2(mu[0], mu[1]))
    report = dict(
        trial=T.TRIAL, source=os.path.basename(T.GRF_RAW_MOT),
        rotation_about_z_deg=float(np.degrees(theta)),
        mean_grf=[float(v) for v in mu],
        mean_grf_magnitude=float(np.linalg.norm(mu)),
        note="expected treadmill-frame tilt; model gravity is pre-rotated; do NOT de-rotate",
    )
    with open(os.path.join(T.HERE, "grf_frame_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"{T.TRIAL}: GRF frame tilt about z = {np.degrees(theta):+.2f} deg, "
          f"|mean| = {np.linalg.norm(mu):.1f} N (expected; model gravity pre-rotated)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())