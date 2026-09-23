"""MocoInverse on the SQR level-walking data (OpenSim Moco, Python API).

Solves the *muscle redundancy* problem over a long time window by splitting it
into short, overlapping segments -- because a single MocoInverse problem over
20 s at the required 2 ms mesh needs ~29 GB of RAM and CasADi dies with
"bad allocation".

Pipeline
--------
    1. Build a ``ModelProcessor``: load ``SQR_simbody.osim``, attach the
       measured GRF as ``ExternalForce`` objects, convert the Millard2012
       muscles to DeGrooteFregly2016 muscles, and add residual/reserve
       ``CoordinateActuator``s so that every degree of freedom is actuated.
    2. Disable the model's own ``HuntCrossleyForce`` foot contacts (the measured
       GRF is prescribed instead, so they would double-count the ground force).
    3. Split ``[--start, --end]`` into slices of ``--segment`` seconds.  For
       each slice solve ``[t_k - overlap, t_{k+1} + overlap]`` (clamped to the
       available data) but KEEP only ``[t_k, t_{k+1}]``.  The extra ``overlap``
       on both sides is *burn-in*: muscle activation has memory
       (tau_act = 10 ms, tau_deact = 40 ms for DeGrooteFregly2016), so the first
       ~5*tau of a solve is contaminated by the unknown initial activation.
       Discarding it makes the kept part insensitive to where the timeline was
       cut, and the retained slices still tile [start, end] exactly.
    4. Stitch the retained slices onto one uniform time grid, and write the
       combined solution plus a seam-continuity report that quantifies whatever
       mismatch is left at the cuts.

Run with the ``opensim_scripting`` conda environment::

    python experiments\\moco\\moco_inverse.py --start 180 --end 200
    ... --segment 2 --overlap 0.5     # defaults
    ... --segment 0                   # single-window mode (small windows only)
    ... --smoke                       # 0.3 s window, to verify the stack

Outputs (in ``experiments/moco/outputs/``):
    moco_inverse_solution.sto        stitched trajectory (states + controls)
    moco_inverse_activations.sto     muscle activations only
    moco_inverse_excitations.sto     muscle excitations only
    moco_inverse_residuals.csv       reserve/residual control magnitudes
    moco_inverse_activation_report.csv
    moco_inverse_seams.csv           seam-continuity diagnostics
    moco_inverse_summary.txt         report incl. per-segment timings
    segments/seg_XX.sto              raw solution of every segment

Notes
-----
* The default mesh interval is 2 ms.  MocoInverse's cost only penalises the
  *magnitude* of the excitation, never its rate of change, and the collocation
  of the activation dynamics admits a mesh-scale zigzag whenever the mesh
  interval h is not small compared with tau.  At h = 20 ms (h/tau = 2) the
  muscle activation itself picks up several percent of non-physiological
  ripple; at 2 ms (h/tau = 0.2) it is ~20x smaller.  Run
  ``check_solution_quality.py`` on the result before drawing conclusions.
* --mesh 0.001 does NOT converge with this Moco build (forward-finite-difference
  conditioning limit), so 0.002 is the practical finest mesh.
"""

import argparse
import os
import re
import sys
import time

import numpy as np
from scipy.signal import butter, filtfilt

# Long segmented runs print a lot; keep the log usable even if the process dies.
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

# Locate the CasADi solver plugins before importing opensim, so that the script
# works even when the conda environment was not activated.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moco_bootstrap                                    # noqa: E402
moco_bootstrap.ensure_casadi_plugins()

import opensim as osim                                   # noqa: E402

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "experiments", "data", "SQR_walking")
OUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
SEG_DIR = os.path.join(OUT_DIR, "segments")

MODEL_FILE = os.path.join(DATA_DIR, "SQR_simbody.osim")
# The 15 Hz filtered kinematics (ik_filter.py), NOT level_walking_ik.mot: that
# one is the raw InverseKinematicsTool output.
IK_FILE = os.path.join(DATA_DIR, "level_walking_ik_filtered.mot")
GRF_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")

# Motion that needs an actuator but has no muscle: the 6 pelvis DOF and the
# 3 lumbar DOF (the torso is connected to the pelvis by a CustomJoint).
#
# Running the existing inverse dynamics on this dataset shows the pelvis needs
# up to ~132 N (tx), ~98 N (ty) and ~78 N*m (tilt) during level walking.  The
# OpenSim Moco examples use (250 N*m, 50 N) here and 50 N is too small for this
# subject -- the translational residuals saturate and the problem becomes
# infeasible.  Keep these as tight as is comfortably feasible: an oversized
# residual lets the solver hide the real dynamics inside the residual.
RESIDUAL_ROTATIONAL_FORCE = 400.0     # N*m, for pelvis_tilt / list / rotation
RESIDUAL_TRANSLATIONAL_FORCE = 500.0  # N,   for pelvis_tx / ty / tz

# Kinematics low-pass cutoff.  MocoInverse fits a spline through the
# coordinates and differentiates it twice, so unfiltered IK noise turns into
# enormous joint accelerations.  This MUST equal the GRF cutoff: inverse
# dynamics / MocoInverse combine twice-differentiated coordinates with the
# external loads, so if the two are band-limited differently the residual does
# not cancel.  Both are 15 Hz, produced by experiments/scripts/ik_filter.py
# (KIN_FILTER_FC) and experiments/scripts/grf_filter.py (GRF_FILTER_FC) and
# configured in Setup_InverseDynamics.xml.
# NOTE: IK_FILE is already that filtered file, so the filter in
# prepare_kinematics_radians() is a second, harmless pass; it is kept so the
# pipeline stays correct if a raw IK file is passed with --filter-hz.
KINEMATICS_FILTER_HZ = 15.0

# Mesh intervals above which a single solve refuses to start (memory guard).
# Measured: ~2.9 MB RAM per interval on this machine, so 3000 intervals
# ~= 8.5 GB; 10,000 intervals (20 s at 2 ms) ~= 28 GB -> "bad allocation".
MAX_MESH_INTERVALS = 3000

# Segmentation defaults.  The activation memory is tau_deact = 40 ms, so an
# overlap of 0.5 s = 12.5 tau is ample; the seam report measures what is left.
DEFAULT_SEGMENT = 2.0     # s kept per segment in the stitched result
DEFAULT_OVERLAP = 0.5     # s of burn-in solved and discarded on each side
GRF_MARGIN = 1.0          # s of extra GRF data kept around the solve span

# When comparing two adjacent segments at a seam, ignore this much of each
# segment's own solve window: the start has an initial-activation transient and
# the end has a free-final-state boundary layer (both measured to die out within
# ~0.1-0.2 s, i.e. a few tau).  0.4 s = 10 tau_deact.
SEAM_SAFE_MARGIN = 0.4


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------
def build_model_processor(grf_xml, use_dgf_muscles=True,
                          add_residuals=True, add_reserves=True,
                          residual_rot=RESIDUAL_ROTATIONAL_FORCE,
                          residual_trans=RESIDUAL_TRANSLATIONAL_FORCE):
    """Return a configured ``osim.ModelProcessor`` for the SQR walking model.

    The operator order matters: external loads and muscle conversion happen
    first, residual/reserve actuators last so that they see the final muscle
    set.
    """
    mp = osim.ModelProcessor(MODEL_FILE)

    # Measured ground reaction forces, applied to calcn_r / calcn_l.
    # ModOpAddExternalLoads ignores the ExternalLoads <datafile> -> model
    # kinematics coupling, so the absolute path inside the XML is used as-is.
    mp.append(osim.ModOpAddExternalLoads(grf_xml))

    if use_dgf_muscles:
        # Rigid tendons + no passive fiber forces = the standard, well
        # conditioned Moco muscle model.
        mp.append(osim.ModOpIgnoreTendonCompliance())
        mp.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
        mp.append(osim.ModOpIgnorePassiveFiberForcesDGF())
        mp.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))

    # Residuals: CoordinateActuators on coordinates of joints whose *parent*
    # body is Ground, i.e. the 6 pelvis DOF.
    if add_residuals:
        mp.append(osim.ModOpAddResiduals(residual_rot, residual_trans, 1.0))

    # Reserves: weak CoordinateActuators on every coordinate that does not
    # already have one (includes the 3 lumbar DOF, which no muscle crosses).
    if add_reserves:
        mp.append(osim.ModOpAddReserves(1.0))

    return mp


# ---------------------------------------------------------------------------
# Text-table helpers
# ---------------------------------------------------------------------------
def read_mot(path):
    """Read an OpenSim .mot/.sto table -> (header_lines, label_line, ndarray)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    i = 0
    header = []
    while i < len(lines) and not lines[i].lstrip().lower().startswith("endheader"):
        header.append(lines[i].rstrip("\r\n"))
        i += 1
    header.append("endheader")
    i += 1
    labels = lines[i].rstrip("\r\n")
    i += 1
    data = np.loadtxt([ln for ln in lines[i:] if ln.strip()])
    return header, labels, data


def write_mot(path, header, labels, data):
    """Write an OpenSim .mot table, fixing up the nRows field."""
    out_header = []
    for h in header:
        if h.strip().lower().startswith("nrows"):
            out_header.append(f"nRows={data.shape[0]}")
        else:
            out_header.append(h)
    with open(path, "w", encoding="utf-8") as f:
        for h in out_header:
            f.write(h + "\n")
        f.write(labels + "\n")
        for row in data:
            f.write("\t".join(f"{v:.8f}" for v in row) + "\n")


def data_time_range(xml_or_mot):
    """Time span of a GRF .mot file (or of the file an ExternalLoads XML names).

    Only the first and last data rows are parsed, so this is cheap even for the
    250k-row raw recording.
    """
    path = xml_or_mot
    if path.lower().endswith(".xml"):
        with open(path, "r", encoding="utf-8") as f:
            m = re.search(r"<datafile>\s*(.*?)\s*</datafile>", f.read(), re.S)
        if not m:
            return None
        path = m.group(1).strip().replace("/", os.sep)
    if not os.path.exists(path):
        return None

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    i_end = next((i for i, l in enumerate(lines)
                  if l.strip().lower().startswith("endheader")), None)
    if i_end is None:
        return None
    i_lab = next((k for k in range(i_end + 1, len(lines)) if lines[k].strip()), None)
    if i_lab is None:
        return None
    data_rows = [l for l in lines[i_lab + 1:] if l.strip()]
    if not data_rows:
        return None
    return float(data_rows[0].split()[0]), float(data_rows[-1].split()[0])


def prepare_grf_window(src_xml, start, end, out_dir, margin=GRF_MARGIN):
    """Write a window-limited copy of the GRF data and a matching loads XML.

    The raw recording is 249.6 s / 249,600 rows (~54 MB of text) and Moco
    re-parses it on every model initialization while building the problem
    (>100 full parses = minutes of pure file I/O).  Trimming it to the span we
    actually solve does not change the solution at all.
    """
    with open(src_xml, "r", encoding="utf-8") as f:
        xml = f.read()
    m = re.search(r"<datafile>\s*(.*?)\s*</datafile>", xml, re.S)
    if not m:
        print("  (no <datafile> in external loads XML; using it as-is)")
        return src_xml
    src_mot = m.group(1).strip().replace("/", os.sep)
    if not os.path.exists(src_mot):
        print(f"  (GRF data file not found: {src_mot}; using XML as-is)")
        return src_xml

    os.makedirs(out_dir, exist_ok=True)
    header, labels, data = read_mot(src_mot)
    t = data[:, 0]
    sel = (t >= start - margin) & (t <= end + margin)
    if not sel.any():
        print(f"  (no GRF samples in {start}..{end}; using XML as-is)")
        return src_xml
    sub = data[sel]

    tag = f"{start:g}_{end:g}".replace(".", "p")
    out_mot = os.path.join(out_dir, f"grf_window_{tag}.mot")
    out_xml = os.path.join(out_dir, f"grf_window_{tag}.xml")
    write_mot(out_mot, header, labels, sub)

    xml_new = re.sub(r"(<datafile>\s*).*?(\s*</datafile>)",
                     lambda mm: mm.group(1) + out_mot.replace(os.sep, "/") + mm.group(2),
                     xml, flags=re.S)
    with open(out_xml, "w", encoding="utf-8") as f:
        f.write(xml_new)

    print(f"  GRF trimmed: {data.shape[0]} -> {sub.shape[0]} rows "
          f"(t {sub[0, 0]:.3f}..{sub[-1, 0]:.3f} s)  -> {os.path.basename(out_mot)}")
    return out_xml


# ---------------------------------------------------------------------------
# Degrees -> radians conversion for the kinematics
# ---------------------------------------------------------------------------
# ``level_walking_ik_filtered.mot`` carries ``inDegrees=yes``: every rotational
# coordinate is stored in degrees (e.g. knee_extension_r = -67.0).
# ``osim.TimeSeriesTable`` passes those numbers through *verbatim* -- it does
# not honour the inDegrees flag -- and Moco interprets coordinate values in
# SimTK internal units, i.e. radians.  Feeding the file in unchanged therefore
# asks the solver for a 67 radian knee angle, ~57x too large, and the problem
# cannot converge.  Moco's own examples avoid this by shipping a radian copy of
# the coordinates; we generate one here.
def prepare_kinematics_radians(model, src_ik, out_dir, filter_hz=KINEMATICS_FILTER_HZ,
                               tag="ik_radians"):
    """Write a filtered, radians-only copy of the IK file (``inDegrees=no``).

    Rotational coordinates are identified from *the model*, not from their
    names, so ``pelvis_tx/ty/tz`` are correctly left alone.
    ``filter_hz <= 0`` disables the low-pass filter.
    """
    header, labels, data = read_mot(src_ik)
    col_names = labels.split()

    rot_coords = set()
    for i in range(model.getCoordinateSet().getSize()):
        coord = model.getCoordinateSet().get(i)
        if coord.getMotionType() == osim.Coordinate.Rotational:
            rot_coords.add(coord.getName())

    t = data[:, 0]
    fs = 1.0 / float(np.median(np.diff(t)))

    if filter_hz and filter_hz > 0:
        if filter_hz >= fs / 2.0:
            print(f"  kinematics: cutoff {filter_hz} Hz >= Nyquist {fs/2:.0f} Hz, "
                  f"skipping filter")
        else:
            b, a = butter(4, filter_hz / (fs / 2.0), btype="low")
            for j in range(1, data.shape[1]):
                data[:, j] = filtfilt(b, a, data[:, j])
            print(f"  kinematics: zero-phase 4th-order Butterworth low-pass at "
                  f"{filter_hz} Hz (fs = {fs:.0f} Hz)")

    converted = []
    for j, name in enumerate(col_names):
        if j == 0:                     # time column
            continue
        if name in rot_coords:
            data[:, j] = np.deg2rad(data[:, j])
            converted.append(name)

    header = [re.sub(r"(?i)inDegrees\s*=\s*yes", "inDegrees=no", h) for h in header]
    out_ik = os.path.join(out_dir, f"{tag}.mot")
    write_mot(out_ik, header, labels, data)

    print(f"  kinematics: converted {len(converted)} rotational coordinate(s) "
          f"deg -> rad (inDegrees=no)  -> {os.path.basename(out_ik)}")
    return out_ik


def disable_contact_forces(model):
    """Turn off the model's HuntCrossleyForce foot contacts, in place."""
    disabled = []
    forces = model.updForceSet()
    for i in range(forces.getSize()):
        force = forces.get(i)
        if "HuntCrossley" in force.getConcreteClassName():
            force.set_appliesForce(False)
            disabled.append(force.getName())
    return disabled


# ---------------------------------------------------------------------------
# Segment planning / solving / stitching
# ---------------------------------------------------------------------------
def build_segment_plan(t0, t1, seg_len, overlap, data_start, data_end):
    """Split ``[t0, t1]`` into retained slices with a burn-in solve window.

    Returns a list of dicts with

    ``retain = (a, b)``  the part kept in the stitched result; the retained
                         slices tile ``[t0, t1]`` exactly, without gaps;
    ``solve  = (s, e)``  what is handed to MocoInverse: the retained slice
                         extended by ``overlap`` on both sides, clamped to the
                         available data range.

    A single-slice plan is solved on ``[t0, t1]`` exactly -- there are no seams,
    so no burn-in is needed.
    """
    if seg_len <= 0 or (t1 - t0) <= seg_len + 1e-9:
        return [dict(k=0, retain=(t0, t1), solve=(t0, t1))]

    plan, a, k = [], t0, 0
    while a < t1 - 1e-9:
        b = min(a + seg_len, t1)
        plan.append(dict(k=k, retain=(a, b),
                         solve=(max(a - overlap, data_start),
                                min(b + overlap, data_end))))
        a, k = b, k + 1
    return plan


def traj_times(sol):
    """Time vector of a MocoTrajectory (SWIG Vector is not numpy-convertible)."""
    n = sol.getNumTimes()
    return np.array([sol.getTime()[i] for i in range(n)])


def export_trajectory(sol):
    """-> (times, state_labels, states, control_labels, controls) as arrays."""
    st = sol.exportToStatesTable()
    ct = sol.exportToControlsTable()
    st_labels = list(st.getColumnLabels())
    ct_labels = list(ct.getColumnLabels())
    states = np.column_stack([st.getDependentColumn(l).to_numpy() for l in st_labels])
    ctrls = np.column_stack([ct.getDependentColumn(l).to_numpy() for l in ct_labels])
    return traj_times(sol), st_labels, states, ct_labels, ctrls


def _interp(grid, ts, data):
    """Linear interpolation of a (nT, nCol) array onto `grid`, column by column."""
    out = np.empty((len(grid), data.shape[1]))
    for j in range(data.shape[1]):
        out[:, j] = np.interp(grid, ts, data[:, j])
    return out


def owners_for(grid, plan):
    """Index of the segment that owns each grid point (retained slices tile)."""
    cuts = np.array([p["retain"][0] for p in plan] + [plan[-1]["retain"][1]])
    return np.clip(np.searchsorted(cuts, grid, side="right") - 1, 0, len(plan) - 1)


def stitch(data, plan, grid):
    """Concatenate retained slices; measure the mismatch left at each seam."""
    _, st_labels, _, ct_labels, _ = data[0]
    states = np.full((len(grid), len(st_labels)), np.nan)
    ctrls = np.full((len(grid), len(ct_labels)), np.nan)
    owner = owners_for(grid, plan)
    for k, (ts, _, S, _, C) in enumerate(data):
        m = owner == k
        if not m.any():
            continue
        states[m] = _interp(grid[m], ts, S)
        ctrls[m] = _interp(grid[m], ts, C)

    act_idx = [i for i, l in enumerate(st_labels) if l.endswith("/activation")]
    exc_idx = [i for i, l in enumerate(ct_labels)
               if not l.rsplit("/", 1)[-1].startswith(("residual_", "reserve_"))]

    # Between two adjacent segments the solve windows overlap; the stitched
    # result takes each point from the segment that owns it.  To measure the
    # stitching error we must compare the two only where BOTH are trustworthy:
    #   * a segment's start carries an initial-activation transient
    #     (MocoInverse imposes a(0) = x(0)); measured decay: 0.96 -> 0.005
    #     within 0.1 s, ~1e-4 after 0.2 s (tau_act = 10 ms);
    #   * a segment's END carries a boundary layer from the free-final-state
    #     transversality condition (measured: up to 0.3 in the last ~0.1 s).
    # So use a window symmetric about the cut whose half-width leaves SAFE
    # seconds of margin to each segment's own solve boundary.
    seams = []
    for k in range(1, len(plan)):
        cut = plan[k]["retain"][0]
        hist = cut - plan[k]["solve"][0]           # history segment k has
        ahead = plan[k - 1]["solve"][1] - cut      # future segment k-1 has
        hw = min(hist, ahead) - SEAM_SAFE_MARGIN
        if hw < 0:
            continue
        m = (grid >= cut - hw) & (grid <= cut + hw)
        if m.sum() < 2:
            continue
        g = grid[m]
        tsA, _, SA, _, CA = data[k - 1]
        tsB, _, SB, _, CB = data[k]
        dS = _interp(g, tsB, SB) - _interp(g, tsA, SA)
        dC = _interp(g, tsB, CB) - _interp(g, tsA, CA)
        i_cut = int(np.argmin(np.abs(grid - cut)))
        jump = (np.abs(states[i_cut, act_idx] - states[i_cut - 1, act_idx])
                if i_cut > 0 else np.zeros(len(act_idx)))
        seams.append(dict(
            cut=cut, half_width=hw, n=int(m.sum()),
            act_rms=float(np.sqrt(np.mean(dS[:, act_idx] ** 2))),
            act_max=float(np.max(np.abs(dS[:, act_idx]))),
            exc_rms=float(np.sqrt(np.mean(dC[:, exc_idx] ** 2))),
            exc_max=float(np.max(np.abs(dC[:, exc_idx]))),
            jump_max=float(np.max(jump)),
        ))
    return st_labels, states, ct_labels, ctrls, seams, act_idx, exc_idx


def _header_keyvals(template):
    d = {}
    for h in template:
        if "=" in h:
            k, v = h.split("=", 1)
            d[k.strip().lower()] = v.strip()
    return d


def write_stitched_sto(path, template, grid, labels, columns, header_override=None,
                       n_states=None, n_controls=None):
    """Write a stitched .sto, reusing a real Moco-written header template.

    Reusing the header of a Moco-written segment file keeps the file readable
    both by ``osim.MocoTrajectory`` and by the pure-MATLAB reader.

    ``n_states`` / ``n_controls`` override the counts that come with the
    template.  Moco requires
    ``num_states + num_controls + num_input_controls + num_multipliers +
    num_derivatives + num_slacks + num_parameters == number of columns``; a file
    whose header disagrees with its columns is unreadable by
    ``osim.MocoTrajectory`` (this bit the activation-only / excitation-only
    exports, which inherited the full-trajectory counts).  The check below
    makes that failure loud instead of silent.
    """
    override = dict(header_override or {})
    if n_states is not None:
        override["num_states"] = n_states
    if n_controls is not None:
        override["num_controls"] = n_controls

    tmpl_vals = _header_keyvals(template)
    extra = sum(int(override.get(k, tmpl_vals.get(k, 0)))
                for k in ("num_input_controls", "num_multipliers",
                          "num_derivatives", "num_slacks", "num_parameters"))
    ns = int(override.get("num_states", tmpl_vals.get("num_states", 0)))
    nc = int(override.get("num_controls", tmpl_vals.get("num_controls", 0)))
    if ns + nc + extra != len(labels):
        raise ValueError(
            f"{os.path.basename(path)}: header counts do not match the columns "
            f"(num_states={ns} + num_controls={nc} + {extra} extra != "
            f"{len(labels)} columns). Pass n_states/n_controls explicitly.")

    with open(path, "w", encoding="utf-8") as f:
        for h in template:
            key = h.split("=")[0].strip().lower()
            if key in override:
                f.write(f"{h.split('=')[0].strip()}={override[key]}\n")
            else:
                f.write(h + "\n")
        f.write("\t".join(["time"] + list(labels)) + "\n")
        for i in range(len(grid)):
            f.write(f"{grid[i]:.8f}\t"
                    + "\t".join(f"{v:.8f}" for v in columns[i]) + "\n")


def sto_header_template(path):
    """Header lines (up to and including ``endheader``) of a .sto file."""
    out = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            out.append(line.rstrip("\r\n"))
            if line.strip().lower().startswith("endheader"):
                break
    return out


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------
def actuator_inventory(model):
    """Map each actuator name -> its maximal force/torque (control scaling)."""
    inv = {}
    for i in range(model.getActuators().getSize()):
        act = model.getActuators().get(i)
        opt = float("nan")
        try:
            ca = osim.CoordinateActuator.safeDownCast(act)
            if ca is not None:
                opt = float(ca.getOptimalForce())
            else:
                mus = osim.Muscle.safeDownCast(act)
                if mus is not None:
                    opt = float(mus.getMaxIsometricForce())
        except Exception:
            pass
        inv[act.getName()] = opt
    return inv


def summarize_actuators(model):
    """Split the model's actuators into muscles / residuals / reserves."""
    muscle_names = {model.getMuscles().get(i).getName()
                    for i in range(model.getMuscles().getSize())}
    names = [model.getActuators().get(i).getName()
             for i in range(model.getActuators().getSize())]
    muscles = [n for n in names if n in muscle_names]
    residuals = [n for n in names if n.startswith("residual_")]
    reserves = [n for n in names if n.startswith("reserve_")]
    other = [n for n in names if n not in muscles + residuals + reserves]
    return names, muscles, residuals, reserves, other


def report_control_usage(ctrl_labels, ctrls, model, out_csv):
    """Peak control / force of every residual and reserve actuator."""
    inv = actuator_inventory(model)
    rows = []
    for j, name in enumerate(ctrl_labels):
        short = name.rsplit("/", 1)[-1]
        if not (short.startswith("residual_") or short.startswith("reserve_")):
            continue
        peak = float(np.max(np.abs(ctrls[:, j])))
        opt = inv.get(short, float("nan"))
        rows.append((short, peak, opt, peak * opt if opt == opt else float("nan")))

    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("actuator,peak_control,optimal_force,peak_force\n")
        for short, peak, opt, force in rows:
            f.write(f"{short},{peak:.6f},{opt:.6f},{force:.6f}\n")
    return rows


def report_muscle_activation(st_labels, states, out_csv=None):
    """Peak and mean activation of every muscle."""
    rows = []
    for j, label in enumerate(st_labels):
        if not label.endswith("/activation"):
            continue
        name = label[len("/forceset/"):-len("/activation")] \
            if label.startswith("/forceset/") else label
        v = states[:, j]
        rows.append((name, float(v.max()), float(v.mean())))
    rows.sort(key=lambda r: -r[1])
    n_sat = sum(1 for r in rows if r[1] > 0.98)
    if out_csv:
        with open(out_csv, "w", encoding="utf-8") as f:
            f.write("muscle,peak_activation,mean_activation\n")
            for name, mx, mean in rows:
                f.write(f"{name},{mx:.6f},{mean:.6f}\n")
    return rows, n_sat


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run MocoInverse on the SQR level-walking data, segment by "
                    "segment, and stitch the result into one trajectory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--start", type=float, default=180.0,
                   help="initial time [s] (data span is 180-240 s)")
    p.add_argument("--end", type=float, default=200.0,
                   help="final time [s]")
    p.add_argument("--segment", type=float, default=DEFAULT_SEGMENT,
                   help="seconds of result kept per segment; 0 disables "
                        "segmentation (single window -- small windows only)")
    p.add_argument("--overlap", type=float, default=DEFAULT_OVERLAP,
                   help="burn-in [s] solved and discarded on each side of a "
                        "segment; must exceed ~5x the 40 ms deactivation time "
                        "constant, i.e. >= 0.2 s")
    p.add_argument("--mesh", type=float, default=0.002,
                   help="mesh interval [s]. Default 2 ms (5 steps per ~10 ms "
                        "activation time constant). 0.02 is only for quick "
                        "exploration -- it produces non-physiological "
                        "excitation/activation chatter")
    p.add_argument("--full", action="store_true",
                   help="use the whole 60 s window (180-240 s)")
    p.add_argument("--smoke", action="store_true",
                   help="tiny 0.3 s window, to verify the stack")
    p.add_argument("--no-dgf", action="store_true",
                   help="keep the original Millard2012 muscles")
    p.add_argument("--keep-contacts", action="store_true",
                   help="do NOT disable the HuntCrossleyForce foot contacts")

    p.add_argument("--no-trim-grf", action="store_true",
                   help="do not trim the GRF file (much slower)")
    p.add_argument("--filter-hz", type=float, default=KINEMATICS_FILTER_HZ,
                   help="low-pass cutoff [Hz] for the IK coordinates; "
                        "0 disables filtering")
    p.add_argument("--residual-rot", type=float, default=RESIDUAL_ROTATIONAL_FORCE,
                   help="pelvis rotational residual capacity [N*m]")
    p.add_argument("--residual-trans", type=float, default=RESIDUAL_TRANSLATIONAL_FORCE,
                   help="pelvis translational residual capacity [N]")
    p.add_argument("--tol", type=float, default=1e-3,
                   help="solver convergence tolerance")
    p.add_argument("--max-iter", type=int, default=3000,
                   help="maximum solver iterations")
    p.add_argument("--allow-large", action="store_true",
                   help="override the mesh-interval memory guard and try anyway")
    p.add_argument("--no-segment-files", action="store_true",
                   help="do not write the raw per-segment solution files")
    p.add_argument("--tag", default=None,
                   help="suffix for output file names")
    return p.parse_args(argv)


def solve_segment(args, mp, ik_file, seg, n_total):
    """Run MocoInverse on one segment.

    Returns ``(solution, keepalive, wall_seconds)``.  ``keepalive`` must be kept
    alive by the caller: ``MocoInverseSolution.getMocoSolution()`` hands back a
    *reference* into the MocoInverseSolution (and the solution derives from the
    model), so dropping the owner produces a dangling reference and the next
    touch is an access violation (0xC0000005).  A fresh model is processed for
    every segment rather than reusing one Model object across solves.
    """
    t0, t1 = seg["solve"]

    model = mp.process()
    if not args.keep_contacts:
        disable_contact_forces(model)
    model.initSystem()

    inverse = osim.MocoInverse()
    inverse.setName(f"SQR_MocoInverse_seg{seg['k']:02d}")
    inverse.setModel(osim.ModelProcessor(model))
    inverse.setKinematics(osim.TableProcessor(ik_file))
    inverse.set_initial_time(t0)
    inverse.set_final_time(t1)
    inverse.set_mesh_interval(args.mesh)
    inverse.set_kinematics_allow_extra_columns(True)
    inverse.set_convergence_tolerance(args.tol)
    inverse.set_constraint_tolerance(args.tol)
    inverse.set_max_iterations(args.max_iter)

    print(f"\n--- segment {seg['k'] + 1}/{n_total}: solve [{t0:.3f}, {t1:.3f}] s "
          f"({t1 - t0:.2f} s, {int(round((t1 - t0) / args.mesh))} intervals), "
          f"keep [{seg['retain'][0]:.3f}, {seg['retain'][1]:.3f}] s")
    t_start = time.time()
    try:
        result = inverse.solve()
    except RuntimeError as exc:
        if "allocation" in str(exc).lower():
            sys.exit(f"\n!! CasADi ran out of memory on segment {seg['k']}.\n"
                     f"   Reduce --segment or --overlap.")
        raise
    elapsed = time.time() - t_start
    sol = result.getMocoSolution()
    print(f"    solved={sol.success()}  iterations={sol.getNumIterations()}  "
          f"objective={sol.getObjective():.6g}  wall={elapsed:.1f} s")
    return sol, (result, model, inverse), elapsed


def main(argv=None):
    args = parse_args(argv)

    if args.full:
        args.start, args.end = 180.0, 240.0
    if args.smoke:
        # A 0.3 s window at the 2 ms production mesh (150 intervals) is a
        # surprisingly hard problem and IPOPT diverges on it; smoke tests only
        # check that the stack works, so use the coarse mesh here.
        args.start, args.end, args.mesh = 180.0, 180.3, 0.02
    if not (args.start < args.end):
        sys.exit("--start must be smaller than --end")

    grf_xml = GRF_XML
    for path in (MODEL_FILE, IK_FILE, grf_xml):
        if not os.path.exists(path):
            sys.exit(f"Missing input file: {path}")
    os.makedirs(OUT_DIR, exist_ok=True)

    grf_span = data_time_range(grf_xml)
    ik_span = data_time_range(IK_FILE)
    data_start = max(grf_span[0] if grf_span else args.start,
                     ik_span[0] if ik_span else args.start)
    data_end = min(grf_span[1] if grf_span else args.end,
                   ik_span[1] if ik_span else args.end)

    plan = build_segment_plan(args.start, args.end, args.segment, args.overlap,
                              data_start, data_end)

    print("=" * 72)
    print("OpenSim", osim.GetVersion(), "| Moco", osim.GetMocoVersion())
    print("CasADi solver available:", osim.MocoCasADiSolver.isAvailable())
    print("=" * 72)
    print(f"model      : {MODEL_FILE}")
    print(f"kinematics : {IK_FILE}")
    print(f"grf        : {grf_xml}")
    print(f"window     : {args.start} - {args.end} s")
    print(f"muscles    : {'Millard2012 (unchanged)' if args.no_dgf else 'DeGrooteFregly2016'}")
    print(f"contacts   : {'KEPT (double-counts GRF!)' if args.keep_contacts else 'disabled'}")
    print(f"data span  : {data_start:.1f} - {data_end:.1f} s")
    print(f"segments   : {len(plan)} x {args.segment:g} s kept, "
          f"{args.overlap:g} s burn-in each side, mesh {args.mesh:g} s")
    if len(plan) > 1 and args.overlap < 0.2:
        print("  WARNING: --overlap < 0.2 s is below ~5x the deactivation time")
        print("           constant (40 ms); each segment's initial-activation")
        print("           transient will leak into the kept part.")

    # ---- memory guard: based on the largest SEGMENT, not the whole window ---
    solve_lens = [s["solve"][1] - s["solve"][0] for s in plan]
    max_intervals = int(round(max(solve_lens) / args.mesh))
    est_gb = max_intervals * 2.9 / 1024.0
    print(f"    largest segment: {max(solve_lens):.2f} s = {max_intervals} "
          f"intervals, heuristic RAM need ~{est_gb:.1f} GB")
    if max_intervals > MAX_MESH_INTERVALS and not args.allow_large:
        sys.exit(
            f"\n!! a single segment needs {max_intervals} mesh intervals "
            f"(~{est_gb:.0f} GB) and CasADi will crash with \"bad allocation\".\n"
            f"   Options: smaller --segment / --overlap, or coarser --mesh\n"
            f"   (--allow-large to try anyway)")
    if args.no_dgf:
        print("\nWARNING: --no-dgf keeps Millard2012EquilibriumMuscle, which does not")
        print("         expose the fiber-length state information MocoInverse needs.")
        print("         MocoInverse in practice requires DeGrooteFregly2016 muscles.")

    if not args.no_trim_grf:
        grf_xml = prepare_grf_window(grf_xml, min(s["solve"][0] for s in plan),
                                     max(s["solve"][1] for s in plan), OUT_DIR)
    print("=" * 72)

    # ---- model and kinematics (inspection copy; each segment gets its own) --
    mp = build_model_processor(grf_xml, use_dgf_muscles=not args.no_dgf,
                               residual_rot=args.residual_rot,
                               residual_trans=args.residual_trans)
    model_info = mp.process()
    if not args.keep_contacts:
        disabled = disable_contact_forces(model_info)
        print(f"disabled contact forces: {disabled or '(none found)'}")
    model_info.initSystem()

    names, muscles, residuals, reserves, other = summarize_actuators(model_info)
    print(f"\nprocessed model: {model_info.getNumCoordinates()} coordinates, "
          f"{len(muscles)} muscles, {len(residuals)} residuals, "
          f"{len(reserves)} reserves, {len(other)} other actuators")

    ik_file = prepare_kinematics_radians(model_info, IK_FILE, OUT_DIR,
                                         filter_hz=args.filter_hz)

    # ---- solve every segment ----------------------------------------------
    print("\nsolving ...  (the CasADi/IPOPT log follows)\n")
    t_all0 = time.time()
    sols, timings, keepalive = [], [], []
    for seg in plan:
        sol, ka, el = solve_segment(args, mp, ik_file, seg, len(plan))
        sols.append(sol)
        timings.append(el)
        keepalive.append(ka)      # keeps the MocoInverseSolution and model alive
        if not args.no_segment_files and sol.success():
            os.makedirs(SEG_DIR, exist_ok=True)
            p = os.path.join(SEG_DIR, f"seg_{seg['k']:02d}.sto")
            try:
                sol.write(p)
            except Exception as exc:
                print(f"    (could not write {p}: {exc})")

    total_elapsed = time.time() - t_all0
    ok = all(s.success() for s in sols)
    motion = sum(p["retain"][1] - p["retain"][0] for p in plan)
    print(f"\n>>> total wall time: {total_elapsed:.1f} s "
          f"({total_elapsed/60.0:.2f} min) for {len(plan)} segment(s), "
          f"{motion:.1f} s of motion")
    print(f">>> all segments solved: {ok}")

    # ---- stitch ------------------------------------------------------------
    data = [export_trajectory(s) for s in sols]
    dt_out = args.mesh / 2.0            # Moco writes the trajectory at mesh/2
    grid = np.arange(args.start, args.end + dt_out * 0.5, dt_out)
    st_labels, states, ct_labels, ctrls, seams, act_idx, exc_idx = stitch(
        data, plan, grid)
    print(f"\nstitched onto {len(grid)} points at dt = {dt_out:g} s "
          f"({len(act_idx)} activations, {len(exc_idx)} excitations, "
          f"{len(ct_labels) - len(exc_idx)} residual/reserve controls)")

    # A real Moco-written segment file supplies the .sto header template, so the
    # stitched file stays readable by osim.MocoTrajectory and the MATLAB reader.
    tmpl_path = os.path.join(SEG_DIR, "seg_00.sto")
    if not os.path.exists(tmpl_path):
        tmpl_path = os.path.join(OUT_DIR, "_seg_template.sto")
        sols[0].write(tmpl_path)
    template = sto_header_template(tmpl_path)

    tag = f"_{args.tag}" if args.tag else ""
    sol_file = os.path.join(OUT_DIR, f"moco_inverse_solution{tag}.sto")
    act_file = os.path.join(OUT_DIR, f"moco_inverse_activations{tag}.sto")
    exc_file = os.path.join(OUT_DIR, f"moco_inverse_excitations{tag}.sto")
    res_file = os.path.join(OUT_DIR, f"moco_inverse_residuals{tag}.csv")
    rep_file = os.path.join(OUT_DIR, f"moco_inverse_activation_report{tag}.csv")
    seam_file = os.path.join(OUT_DIR, f"moco_inverse_seams{tag}.csv")
    sum_file = os.path.join(OUT_DIR, f"moco_inverse_summary{tag}.txt")

    objective = float(sum(s.getObjective() for s in sols))
    header_override = {
        "num_iterations": sum(int(s.getNumIterations()) for s in sols),
        "objective": f"{objective:.6f}",
        "objective_excitation_effort": f"{objective:.6f}",
        "solver_duration": f"{float(sum(s.getSolverDuration() for s in sols)):.6f}",
        "status": "Solve_Succeeded" if ok else "Failed",
        "success": "true" if ok else "false",
    }
    write_stitched_sto(sol_file, template, grid, st_labels + ct_labels,
                       np.column_stack([states, ctrls]), header_override,
                       n_states=len(st_labels), n_controls=len(ct_labels))
    write_stitched_sto(act_file, template, grid,
                       [st_labels[i] for i in act_idx], states[:, act_idx],
                       n_states=len(act_idx), n_controls=0)
    write_stitched_sto(exc_file, template, grid,
                       [ct_labels[i] for i in exc_idx], ctrls[:, exc_idx],
                       n_states=0, n_controls=len(exc_idx))
    print(f"wrote stitched solution  -> {sol_file}")
    print(f"wrote activations        -> {act_file}")
    print(f"wrote excitations        -> {exc_file}")

    try:
        chk = osim.MocoTrajectory(sol_file)
        print(f"read-back check: OK ({chk.getNumStates()} states, "
              f"{chk.getNumControls()} controls, {chk.getNumTimes()} times)")
    except Exception as exc:
        print(f"read-back check FAILED ({exc}); the stitched .sto may still be "
              f"readable by the MATLAB reader")

    # ---- reports -----------------------------------------------------------
    rows = report_control_usage(ct_labels, ctrls, model_info, res_file)
    act_rows, n_sat = report_muscle_activation(st_labels, states, rep_file)
    print(f"wrote residual report    -> {res_file}")
    print(f"wrote activation report  -> {rep_file}")

    print(f"\npeak muscle activation:")
    print(f"  {'muscle':<24}{'peak':>8}{'mean':>8}")
    for name, mx, mean in act_rows[:8]:
        flag = "   <-- SATURATED" if mx > 0.98 else ""
        print(f"  {name:<24}{mx:>8.3f}{mean:>8.3f}{flag}")
    print(f"  ... {len(act_rows)} muscles total, {n_sat} saturating at activation 1.0")

    rows_sorted = sorted(rows, key=lambda r: -abs(r[3] if r[3] == r[3] else 0))
    print("\nlargest residual/reserve forces:")
    print(f"  {'actuator':<32}{'peak ctrl':>12}{'optimal':>12}{'peak force':>14}")
    for short, peak, opt, force in rows_sorted[:8]:
        print(f"  {short:<32}{peak:>12.4f}{opt:>12.2f}{force:>14.3f}")

    # ---- seam continuity ---------------------------------------------------
    with open(seam_file, "w", encoding="utf-8") as f:
        f.write("cut_time_s,half_width_s,n_points,act_rms_diff,act_max_diff,"
                "exc_rms_diff,exc_max_diff,stitched_act_jump\n")
        for s in seams:
            f.write(f"{s['cut']:.4f},{s['half_width']:.4f},{s['n']},"
                    f"{s['act_rms']:.6f},{s['act_max']:.6f},"
                    f"{s['exc_rms']:.6f},{s['exc_max']:.6f},{s['jump_max']:.6f}\n")

    print("\nseam continuity  (two adjacent segments compared over a window")
    print(f"                  symmetric about each cut, both >= "
          f"{SEAM_SAFE_MARGIN:g} s from their own solve boundary;")
    print("                  muscle activation units, activation range 0-1):")
    if seams:
        print(f"  {'cut [s]':>9}{'+- [s]':>8}{'n':>6}{'act RMS':>10}{'act max':>10}"
              f"{'exc RMS':>10}{'exc max':>10}{'jump':>9}")
        for s in seams:
            print(f"  {s['cut']:>9.3f}{s['half_width']:>8.2f}{s['n']:>6}"
                  f"{s['act_rms']:>10.5f}{s['act_max']:>10.5f}"
                  f"{s['exc_rms']:>10.5f}{s['exc_max']:>10.5f}{s['jump_max']:>9.5f}")
        print(f"  worst activation mismatch over all seams: "
              f"{max(s['act_max'] for s in seams):.5f}")
        print(f"  worst jump in the stitched activation:    "
              f"{max(s['jump_max'] for s in seams):.5f}")
    else:
        print("  (single segment -- no seams)")
    print(f"wrote seam report        -> {seam_file}")

    # ---- summary -----------------------------------------------------------
    with open(sum_file, "w", encoding="utf-8") as f:
        f.write(f"OpenSim {osim.GetVersion()} / Moco {osim.GetMocoVersion()}\n")
        f.write(f"model={MODEL_FILE}\nkinematics={IK_FILE}\ngrf={grf_xml}\n")
        f.write(f"window={args.start}..{args.end} mesh={args.mesh}\n")
        f.write(f"segments={len(plan)} segment_len={args.segment} "
                f"overlap={args.overlap}\n")
        for p, el in zip(plan, timings):
            f.write(f"  seg{p['k']:02d} solve={p['solve'][0]:.3f}..{p['solve'][1]:.3f}"
                    f" keep={p['retain'][0]:.3f}..{p['retain'][1]:.3f}"
                    f" wall_s={el:.1f}\n")
        f.write(f"wall_time_s={total_elapsed:.1f}\n")
        f.write(f"solved={ok}\n")
        f.write(f"objective={objective:.8g}\n")
        f.write(f"num_states={len(st_labels)}\n")
        f.write(f"num_controls={len(ct_labels)}\n")
        f.write(f"solver_duration_s="
                f"{float(sum(s.getSolverDuration() for s in sols)):.1f}\n")
        f.write(f"solver_iterations={sum(int(s.getNumIterations()) for s in sols)}\n")
        if seams:
            f.write("seam_max_activation_mismatch="
                    f"{max(s['act_max'] for s in seams):.6f}\n")
        f.write(f"activation_saturated_muscles={n_sat}\n")
    print(f"wrote summary            -> {sum_file}")

    print("\nNEXT: run  python experiments/moco/check_solution_quality.py "
          f"\"{sol_file}\"")
    print("      to check the mesh-scale chatter before using these muscle results.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
