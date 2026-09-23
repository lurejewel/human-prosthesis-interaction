"""Residual Reduction Algorithm (RRA) on the SQR level-walking trial.

Entry point: run_ik_rra.py, which runs IK (ik.py) -> 15 Hz filter
(ik_filter.py) -> this module.  RRA cannot produce its own kinematics, so this
file is a library rather than a standalone CLI; grf_filter.py and ik_filter.py
must have run first (the entry point does both).

What it does:
  1. prepare()          builds the muscle-free model, the residual/reserve
                        actuators and the RRA setup XMLs, then trims the GRF and
                        the FILTERED IK to the solve window;
  2. run_rra_passes()   n rounds of RRA.  Round k uses the model produced by
                        round k-1; round 1 uses SQR_RRA_model.osim.  The run
                        stops immediately if a round diverges, contains NaN or
                        saturates a reserve actuator.

RRA is fed the filtered kinematics, never the raw IK (see ik_filter.py), and RRA's
own lowpass is DISABLED (RRA_LOWPASS_FC = -1) because that file has already been
low-passed once: running the same cutoff again would square the response (-6 dB
-> -12 dB at 15 Hz) and quietly narrow the band in which the residuals are
computed.

Generated files (nothing to edit by hand), in data/SQR_walking/RRA_setup/:
    SQR_RRA_model.osim                 muscle-free model + the actuator set
    SQR_RRA_Actuators.xml              the same actuator set, in the reference
                                       format (see below)
    SQR_RRA_Tasks.xml                  CMC_TaskSet: which coordinates are tracked
    grf_window.mot / .xml              external loads trimmed to the window
    level_walking_ik_rra_window.mot    filtered IK trimmed to the window
Results (RRA_pass<k>_*.sto/_adjusted.osim) go to results_dir (ResultsRRA/).

The actuator and task files follow the two reference files in that directory,
SQR_RRA_Actuators_REFERENCE.xml and SQR_RRA_Tasks_REFERENCE.xml, which are the
authority for every name and number in them; the generated files are compared
against those references byte for byte.  In short:
    residuals     FX FY FZ   global PointActuators on the pelvis, acting at its
                             mass center, optimal force 5 / 10 / 5 N
                  MX MY MZ   global TorqueActuators, pelvis vs ground,
                             optimal force 5 / 5 / 10 Nm
                             all six: min/max control = -infinity / infinity
    reserves      one CoordinateActuator per non-pelvis coordinate, named
                  exactly like that coordinate, min/max control -2 / 2, with the
                  optimal forces in COORD_ACTUATOR_FORCE below
    task weights  pelvis_tilt 50; pelvis_list/rotation/t* 20; hip_flexion 50;
                  knee_extension 50; ankle_dorsiflexion 50; hip_ab/adduction and
                  hip_rotation 10; lumbar_* 10
No ControlSet/ControlConstraints file is written: OpenSim calls it deprecated and
"generally unnecessary", and the limits above are carried by the actuators
themselves.  The rest of the RRA parameters:
    kp = 100, kv = 20 (= 2*sqrt(kp), critical damping)
    integrator          error tolerance 1e-3, min step 1e-10, max step 1
    lowpass cutoff      -1 (disabled: the desired kinematics are already
                        filtered by ik_filter.py; the 15 Hz that matters is the
                        one shared with the GRF)
    adjust COM          only on the first pass; the RRA reference always does

Library use:
    import rra
    history = rra.run(ik_filtered=".../level_walking_ik_filtered.mot",
                      t0=170.0, t1=240.0, n_passes=5)
"""
import os
import sys

import numpy as np

import opensim as osim

# make the sibling modules importable when this file is run from anywhere
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grf_filter import read_table, write_table, write_loads_xml      # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data", "SQR_walking")

MODEL_IN = os.path.join(DATA_DIR, "SQR_simbody.osim")
GRF_FILT_MOT = os.path.join(DATA_DIR, "level_walking_grf_filtered.mot")
GRF_FILT_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")

RRA_SETUP_DIR = os.path.join(DATA_DIR, "RRA_setup")
RRA_RESULTS = os.path.join(DATA_DIR, "ResultsRRA")

# file names inside setup_dir / results_dir
MODEL_RRA_NAME = "SQR_RRA_model.osim"
ACTUATORS_XML_NAME = "SQR_RRA_Actuators.xml"
TASKS_XML_NAME = "SQR_RRA_Tasks.xml"
GRF_WINDOW_MOT_NAME = "grf_window.mot"
GRF_WINDOW_XML_NAME = "grf_window.xml"
IK_WINDOW_MOT_NAME = "level_walking_ik_rra_window.mot"

# reference files the generated setup is compared against (never written to)
ACTUATORS_XML_REF_NAME = "SQR_RRA_Actuators_REFERENCE.xml"
TASKS_XML_REF_NAME = "SQR_RRA_Tasks_REFERENCE.xml"

# defaults for the standard directories
MODEL_RRA = os.path.join(RRA_SETUP_DIR, MODEL_RRA_NAME)
ACTUATORS_XML = os.path.join(RRA_SETUP_DIR, ACTUATORS_XML_NAME)
TASKS_XML = os.path.join(RRA_SETUP_DIR, TASKS_XML_NAME)
IK_WINDOW_MOT = os.path.join(RRA_SETUP_DIR, IK_WINDOW_MOT_NAME)
GRF_WINDOW_XML = os.path.join(RRA_SETUP_DIR, GRF_WINDOW_XML_NAME)

DEFAULT_PASSES = 5
DEFAULT_TIME_RANGE = (170.0, 240.0)

# RRA's own lowpass on the desired kinematics.  -1 disables it, which is what we
# want: the window file RRA reads is ik_filter.py's output and is already
# low-passed.  (OpenSim's own documentation of this property: "Low-pass cut-off
# frequency for filtering the desired kinematics.  A negative value results in
# no filtering.  The default value is -1.0, so no filtering.")
RRA_LOWPASS_FC = -1.0

# integrator / optimizer settings, from subject01_Setup_RRA.xml
ERROR_TOLERANCE = 1e-3
MIN_DT = 1e-10
MAX_DT = 1.0
# The RK-Merson integrator at error tolerance 1e-3 takes ~18,000 steps per second
# of simulated time on this model (measured), so 70 s needs ~1.24e6 steps.  The
# OpenSim default of 20,000 would abort after ~1 s; keep generous headroom.
MAX_STEPS_PER_SECOND = 30000

# ---------------------------------------------------------------------------
# Actuator set -- every value below is taken from
# RRA_setup/SQR_RRA_Actuators_REFERENCE.xml and must not drift from it.
# ---------------------------------------------------------------------------
# The six pelvis DOF are driven by residuals laid out as in the OpenSim Gait2354
# RRA example: three global forces acting at the pelvis mass center and three
# global torques.  (name, axis, optimal force)
RESIDUAL_FORCES = (("FX", (1.0, 0.0, 0.0), 5.0),
                   ("FY", (0.0, 1.0, 0.0), 10.0),
                   ("FZ", (0.0, 0.0, 1.0), 5.0))
RESIDUAL_TORQUES = (("MX", (1.0, 0.0, 0.0), 5.0),
                    ("MY", (0.0, 1.0, 0.0), 5.0),
                    ("MZ", (0.0, 0.0, 1.0), 10.0))
RESIDUAL_BODY = "pelvis"            # forces act on the pelvis, at its mass center
RESIDUAL_BODY_B = "ground"          # ... the torques act between pelvis and ground
RESIDUAL_CONTROL = float("inf")     # min/max control: unlimited

# The pelvis DOF above are covered by the residuals, so they get no
# CoordinateActuator; every other coordinate gets exactly one, named like the
# coordinate it drives.  optimal force [N*m], min/max control = +/-2.
#
# NOTE (measured, not a hypothetical): build_model() removes the muscles, so
# these actuators have to produce the whole joint moment themselves, and at
# +/-2 the caps sit below this subject's peaks -- knee_extension_l needs 137 Nm
# against a 100 Nm cap (121 Nm already in 170-175 s) and hip_flexion_l needs
# 257 Nm against 200 Nm.  RRA then stops with "Ipopt: Maximum iterations
# exceeded ... Model cannot generate the forces necessary to achieve the target
# acceleration".  The values are kept as they are because the reference files
# are the authority; widen COORD_ACTUATOR_CONTROL to 3.0 (verified working, no
# other value touched) or keep the muscles in the model if that has to change.
COORD_ACTUATOR_FORCE = {
    "hip_flexion_r": 100.0, "hip_adduction_r": 100.0, "hip_rotation_r": 20.0,
    "knee_extension_r": 50.0, "ankle_dorsiflexion_r": 200.0,
    "hip_flexion_l": 100.0, "hip_adduction_l": 100.0, "hip_rotation_l": 20.0,
    "knee_extension_l": 50.0, "ankle_dorsiflexion_l": 200.0,
    "lumbar_extension": 60.0, "lumbar_bending": 60.0, "lumbar_rotation": 20.0,
}
COORD_ACTUATOR_CONTROL = 2.0
PELVIS_DOF = ("pelvis_tilt", "pelvis_list", "pelvis_rotation",
              "pelvis_tx", "pelvis_ty", "pelvis_tz")

# The control limits live on the actuators themselves, which is why no
# ControlConstraints/ControlSet file is generated any more: OpenSim itself calls
# that file deprecated ("generally unnecessary" in the RRA log) and it only
# repeated what min_control/max_control on the actuators already say.
RESIDUAL_FORCE_NAMES = tuple(n for n, _, _ in RESIDUAL_FORCES)
RESIDUAL_TORQUE_NAMES = tuple(n for n, _, _ in RESIDUAL_TORQUES)
ACTUATOR_NAMES = RESIDUAL_FORCE_NAMES + RESIDUAL_TORQUE_NAMES \
    + tuple(COORD_ACTUATOR_FORCE)

# RRA task weights, from RRA_setup/SQR_RRA_Tasks_REFERENCE.xml
TASK_WEIGHTS = {
    "pelvis_tilt": 50.0,
    "pelvis_list": 20.0, "pelvis_rotation": 20.0,
    "pelvis_tx": 20.0, "pelvis_ty": 20.0, "pelvis_tz": 20.0,
    "hip_flexion_r": 50.0, "hip_adduction_r": 10.0, "hip_rotation_r": 10.0,
    "knee_extension_r": 50.0, "ankle_dorsiflexion_r": 50.0,
    "hip_flexion_l": 50.0, "hip_adduction_l": 10.0, "hip_rotation_l": 10.0,
    "knee_extension_l": 50.0, "ankle_dorsiflexion_l": 50.0,
    "lumbar_extension": 10.0, "lumbar_bending": 10.0, "lumbar_rotation": 10.0,
}
TASK_WEIGHT_DEFAULT = 10.0
KP, KV, KA = 100.0, 20.0, 1.0   # kv = 2*sqrt(kp) -> critically damped

# failure thresholds
RESIDUAL_FORCE_WARN = 25.0      # N
RESIDUAL_TORQUE_WARN = 50.0     # Nm
SATURATION_FRACTION = 0.999     # |control| above this counts as saturated
SATURATION_ALLOWED = 0.05       # fraction of samples allowed saturated
DIVERGENCE_DEG = 720.0          # any coordinate beyond this is divergence
ADJUSTED_COM_BODY = "torso"     # heaviest segment; the RRA reference uses torso


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def read_sto(path):
    """Read an OpenSim .sto whose header may span several lines.

    `read_table` assumes a single label line; RRA's trajectory files put the
    column names after a blank line and may add a 'der' suffix line, so this
    variant locates 'endheader' explicitly.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    hi = next(i for i, l in enumerate(lines) if l.strip() == "endheader")
    j = hi + 1
    while j < len(lines) and not lines[j].strip():
        j += 1
    names = lines[j].split()
    rows = [l for l in lines[j + 1:] if l.strip()]
    return names, np.loadtxt(rows, ndmin=2)


def paths(setup_dir=RRA_SETUP_DIR, results_dir=RRA_RESULTS):
    """The files a prepared RRA run reads and writes, for the given directories.

    Everything RRA needs is either produced by build_model()/write_window_files()
    inside `setup_dir`, or written by RRA itself into `results_dir`.
    """
    return {
        "model": os.path.join(setup_dir, MODEL_RRA_NAME),
        "actuators_xml": os.path.join(setup_dir, ACTUATORS_XML_NAME),
        "tasks_xml": os.path.join(setup_dir, TASKS_XML_NAME),
        "grf_window_mot": os.path.join(setup_dir, GRF_WINDOW_MOT_NAME),
        "grf_window_xml": os.path.join(setup_dir, GRF_WINDOW_XML_NAME),
        "ik_window": os.path.join(setup_dir, IK_WINDOW_MOT_NAME),
        "results_dir": results_dir,
    }


def _num(x):
    """Number exactly as OpenSim writes it in XML: 17 significant digits.  That
    gives '5' / '-2' / '0' for the integral values used here, and reproduces the
    reference's '-0.084115834136295586' for the pelvis mass center (Python's
    repr would shorten it to a 16-digit form that no longer matches)."""
    return "{:.17g}".format(float(x))


def _vec3(v):
    """A SimTK Vec3 as a 3-tuple of floats (SWIG's Vec3 does not iterate)."""
    return (float(v[0]), float(v[1]), float(v[2]))


def add_actuators(model, verbose=True):
    """Add the residual + reserve actuator set described by the constants above.

    This replaces the ModOpAddResiduals()/ModOpAddReserves() route, which named
    the actuators 'residual_...'/'reserve_...' and could not express the
    reference layout (FX/FY/FZ as global PointActuators, MX/MY/MZ as global
    TorqueActuators, one CoordinateActuator per remaining coordinate).

    The residual forces act at the pelvis mass center, so the point follows the
    model: the reference's '-0.084115834136295586 0 0' is that mass center in
    SQR_simbody.osim and stays correct if the model is rescaled.

    Returns the list of actuator names, in file order.
    """
    point = model.getBodySet().get(RESIDUAL_BODY).getMassCenter()

    for name, axis, force in RESIDUAL_FORCES:
        a = osim.PointActuator()
        a.setName(name)
        a.set_body(RESIDUAL_BODY)
        a.set_point(point)
        a.set_point_is_global(False)
        a.set_direction(osim.Vec3(*axis))
        a.set_force_is_global(True)
        a.set_optimal_force(force)
        a.set_min_control(-RESIDUAL_CONTROL)
        a.set_max_control(RESIDUAL_CONTROL)
        model.addForce(a)

    for name, axis, force in RESIDUAL_TORQUES:
        a = osim.TorqueActuator()
        a.setName(name)
        a.set_bodyA(RESIDUAL_BODY)
        a.set_bodyB(RESIDUAL_BODY_B)
        a.set_torque_is_global(True)
        a.set_axis(osim.Vec3(*axis))
        a.set_optimal_force(force)
        a.set_min_control(-RESIDUAL_CONTROL)
        a.set_max_control(RESIDUAL_CONTROL)
        model.addForce(a)

    coords = [model.getCoordinateSet().get(i).getName()
              for i in range(model.getCoordinateSet().getSize())]
    unknown = [c for c in coords
               if c not in PELVIS_DOF and c not in COORD_ACTUATOR_FORCE]
    if unknown:
        sys.exit("no actuator is defined for coordinate(s) {} -- add them to "
                 "COORD_ACTUATOR_FORCE (see SQR_RRA_Actuators_REFERENCE.xml)"
                 .format(", ".join(unknown)))
    for c in coords:
        if c in PELVIS_DOF:
            continue
        a = osim.CoordinateActuator()
        a.setName(c)                        # named exactly like its coordinate
        a.set_coordinate(c)
        a.set_optimal_force(COORD_ACTUATOR_FORCE[c])
        a.set_min_control(-COORD_ACTUATOR_CONTROL)
        a.set_max_control(COORD_ACTUATOR_CONTROL)
        model.addForce(a)

    if verbose:
        print("   actuators: {} residuals (FX FY FZ / MX MY MZ) + {} coordinate "
              "actuators = {}".format(len(RESIDUAL_FORCES) + len(RESIDUAL_TORQUES),
                                      len(COORD_ACTUATOR_FORCE), len(ACTUATOR_NAMES)))
    return list(ACTUATOR_NAMES)


# ---------------------------------------------------------------------------
# Step 1: build the RRA model and setup files
# ---------------------------------------------------------------------------
def build_model(setup_dir=RRA_SETUP_DIR, model_in=MODEL_IN, verbose=True):
    """Muscle-free model + the actuator set above; written to disk.

    RRA integrates the model forward and drives it with a tracking controller,
    so every degree of freedom needs an Actuator.  Muscles are not actuators in
    that sense, hence ModOpRemoveMuscles().

    Returns the path of the generated model.
    """
    os.makedirs(setup_dir, exist_ok=True)
    model_out = os.path.join(setup_dir, MODEL_RRA_NAME)

    mp = osim.ModelProcessor(model_in)
    mp.append(osim.ModOpRemoveMuscles())
    model = mp.process()
    model.setName("SQR_RRA")

    # the model carries an Umberger2010 metabolics probe referring to the muscles
    # we just removed; without this it warns once per muscle on every load
    model.updProbeSet().clearAndDestroy()

    add_actuators(model, verbose=verbose)
    model.initSystem()

    # the model must carry exactly the reference actuator set, in file order.
    # Checked on getActuators(), not on the ForceSet: the base model also brings
    # four passive forces (foot_r/foot_l/knee_lim_r/knee_lim_l) that are not
    # actuators and do not appear in the reference either.
    acts = model.getActuators()
    got = [acts.get(i).getName() for i in range(acts.getSize())]
    if got != list(ACTUATOR_NAMES):
        sys.exit("actuator set/order differs from the reference:\n  got  {}\n  want {}"
                 .format(got, list(ACTUATOR_NAMES)))
    if verbose:
        print("   model: {} bodies, total mass {:.1f} kg, {} forces of which {} "
              "actuators".format(model.getBodySet().getSize(),
                                 model.getTotalMass(model.initSystem()),
                                 model.getForceSet().getSize(), len(got)))
    model.printToXML(model_out)
    if verbose:
        print("   wrote {}".format(model_out))

    write_actuators_xml(model, os.path.join(setup_dir, ACTUATORS_XML_NAME))
    write_tasks_xml(model, os.path.join(setup_dir, TASKS_XML_NAME))
    return model_out


# ---------------------------------------------------------------------------
# SQR_RRA_Actuators.xml -- must match SQR_RRA_Actuators_REFERENCE.xml byte for
# byte, so the indentation below is transcribed from the reference as it is.
# ---------------------------------------------------------------------------
# The reference's three PointActuator blocks do NOT share one indentation: they
# carry the mixed spaces/tabs left over from hand-editing the Gait2354 example.
# (open tag, the body/point_is_global/force_is_global lines, close tag)
_POINT_INDENT = {
    "FX": ("       \t\t", "          \t\t", "        \t"),
    "FY": ("    \t\t", "      \t\t\t", "    \t\t"),
    "FZ": ("      \t\t", "          \t\t", "        \t"),
}
_TAG_INDENT = "\t\t\t"
_TORQUE_BODY_INDENT = "            \t"
_PROP_INDENT = "\t\t\t\t"


def _point_actuator_lines(name, axis, force, point):
    ind, ind2, close = _POINT_INDENT[name]
    return [
        '{}<PointActuator name="{}">'.format(ind, name),
        '{}<isDisabled> false </isDisabled>'.format(_PROP_INDENT),
        '{}<min_control> -infinity </min_control>'.format(_PROP_INDENT),
        '{}<max_control> infinity </max_control>'.format(_PROP_INDENT),
        '{}<body> {} </body>'.format(ind2, RESIDUAL_BODY),
        '{}<point>{}</point>'.format(_PROP_INDENT, point),
        '{}<point_is_global> false </point_is_global>'.format(ind2),
        '{}<direction> {} </direction>'.format(
            _PROP_INDENT, " ".join(_num(v) for v in axis)),
        '{}<force_is_global> true </force_is_global>'.format(ind2),
        '{}<optimal_force> {} </optimal_force>'.format(_PROP_INDENT, _num(force)),
        '{}</PointActuator>'.format(close),
    ]


def _torque_actuator_lines(name, axis, force):
    ind, ind2 = _TAG_INDENT, _TORQUE_BODY_INDENT
    return [
        '{}<TorqueActuator name="{}">'.format(ind, name),
        '{}<isDisabled> false </isDisabled>'.format(_PROP_INDENT),
        '{}<min_control> -infinity </min_control>'.format(_PROP_INDENT),
        '{}<max_control> infinity </max_control>'.format(_PROP_INDENT),
        '{}<bodyA> {} </bodyA>'.format(ind2, RESIDUAL_BODY),
        '{}<bodyB> {} </bodyB>'.format(ind2, RESIDUAL_BODY_B),
        '{}<torque_is_global> true </torque_is_global>'.format(ind2),
        '{}<axis> {} </axis>'.format(
            _PROP_INDENT, " ".join(_num(v) for v in axis)),
        '{}<optimal_force> {} </optimal_force>'.format(_PROP_INDENT, _num(force)),
        '{}</TorqueActuator>'.format(ind),
    ]


def _coord_actuator_lines(name, force):
    ind = _TAG_INDENT
    return [
        '{}<CoordinateActuator name="{}">'.format(ind, name),
        '{}<min_control> {} </min_control>'.format(
            _PROP_INDENT, _num(-COORD_ACTUATOR_CONTROL)),
        '{}<max_control> {} </max_control>'.format(
            _PROP_INDENT, _num(COORD_ACTUATOR_CONTROL)),
        '{}<coordinate>{}</coordinate>'.format(_PROP_INDENT, name),
        '{}<optimal_force>{}</optimal_force>'.format(_PROP_INDENT, _num(force)),
        '{}</CoordinateActuator>'.format(ind),
    ]


def write_actuators_xml(model, path):
    """SQR_RRA_Actuators.xml: the model's actuator set, in the reference format.

    Written from the constants above -- the same ones add_actuators() puts into
    the model -- rather than from OpenSim's printToXML(), which adds the schema
    documentation comments that the reference file does not have.
    """
    point = " ".join(_num(v) for v in _vec3(
        model.getBodySet().get(RESIDUAL_BODY).getMassCenter()))
    lines = ['<?xml version="1.0" encoding="UTF-8" ?>',
             '<OpenSimDocument Version="40600">',
             '\t<ForceSet name="SQR_RRA">',
             '\t\t<objects>']
    for name, axis, force in RESIDUAL_FORCES:
        lines += _point_actuator_lines(name, axis, force, point)
    for name, axis, force in RESIDUAL_TORQUES:
        lines += _torque_actuator_lines(name, axis, force)
    for name, force in COORD_ACTUATOR_FORCE.items():
        lines += _coord_actuator_lines(name, force)
    lines += ['\t\t</objects>', '\t\t<groups />', '\t</ForceSet>',
              '</OpenSimDocument>', '']

    with open(path, "w", encoding="utf-8", newline="\r\n") as f:
        f.write("\n".join(lines))
    print("   wrote {} ({} actuators)".format(path, len(ACTUATOR_NAMES)))


def write_tasks_xml(model, path):
    """SQR_RRA_Tasks.xml: one tracking task per coordinate, as in the reference.

    Format and weights follow RRA_setup/SQR_RRA_Tasks_REFERENCE.xml exactly
    (2-space indentation here, unlike the actuator file).
    """
    coords = [model.getCoordinateSet().get(i).getName()
              for i in range(model.getCoordinateSet().getSize())]
    docs = ['<?xml version="1.0" encoding="UTF-8"?>',
            '<OpenSimDocument Version="40500">',
            '  <CMC_TaskSet name="SQR_RRA">',
            '    <objects>']
    for c in coords:
        w = TASK_WEIGHTS.get(c, TASK_WEIGHT_DEFAULT)
        docs += [
            '      <CMC_Joint name="{}">'.format(c),
            '        <on>true</on>',
            '        <active>true false false</active>',
            '        <weight>{}</weight>'.format(w),
            '        <wrt_body>-1</wrt_body>',
            '        <express_body>-1</express_body>',
            '        <kp>{} 1 1</kp>'.format(KP),
            '        <kv>{} 1 1</kv>'.format(KV),
            '        <ka>{} 1 1</ka>'.format(KA),
            '        <r0>0 0 0</r0>',
            '        <r1>0 0 0</r1>',
            '        <r2>0 0 0</r2>',
            '        <coordinate>{}</coordinate>'.format(c),
            '        <limit>0</limit>',
            '      </CMC_Joint>',
        ]
    docs += ['    </objects>', '    <groups />', '  </CMC_TaskSet>', '</OpenSimDocument>', '']
    with open(path, "w", encoding="utf-8", newline="\r\n") as f:
        f.write("\n".join(docs))
    print("   wrote {} ({} tracking tasks)".format(path, len(coords)))


# ---------------------------------------------------------------------------
# Step 2: window files
# ---------------------------------------------------------------------------
def write_window_files(ik_filtered, t0, t1, setup_dir=RRA_SETUP_DIR, margin=1.0):
    """Trim the GRF and the FILTERED IK to the solve window so RRA does not
    re-parse the 249,600-row GRF on every model initialization.

    `ik_filtered` is the 15 Hz kinematics (ik_filter.py), never the raw IK.
    Returns (grf_window_mot, grf_window_xml, ik_window_mot).
    """
    print("\n   preparing window files for {:.0f}-{:.0f} s ...".format(t0, t1))
    os.makedirs(setup_dir, exist_ok=True)
    grf_win = os.path.join(setup_dir, GRF_WINDOW_MOT_NAME)
    grf_xml = os.path.join(setup_dir, GRF_WINDOW_XML_NAME)
    ik_win = os.path.join(setup_dir, IK_WINDOW_MOT_NAME)

    header, labels, data = read_table(GRF_FILT_MOT)
    sel = (data[:, 0] >= t0 - margin) & (data[:, 0] <= t1 + margin)
    write_table(grf_win, header, labels, data[sel])
    write_loads_xml(grf_xml, grf_win, template_xml=GRF_FILT_XML)
    print("   GRF: {} -> {} rows".format(data.shape[0], int(sel.sum())))

    header, labels, data = read_table(ik_filtered)
    # both bounds are on the time column (column 0); column 1 is pelvis_tilt
    sel = (data[:, 0] >= t0) & (data[:, 0] <= t1)
    if not sel.any():
        sys.exit("no IK samples inside {:.0f}-{:.0f} s".format(t0, t1))
    write_table(ik_win, header, labels, data[sel])
    print("   IK : {} -> {} rows".format(data.shape[0], int(sel.sum())))
    return grf_win, grf_xml, ik_win


# ---------------------------------------------------------------------------
# Step 3: RRA
# ---------------------------------------------------------------------------
def analyse_actuation(path, t0, t1):
    """Peak residual force/torque and the reserve usage, from an Actuation_force.sto.

    The columns are named after the actuators (FX FY FZ / MX MY MZ / coordinate
    names) and hold generalized forces in N or N*m, i.e. control * optimal_force.
    A coordinate actuator may therefore reach optimal_force * COORD_ACTUATOR_CONTROL
    (+/-2), which is the limit the saturation test uses.
    """
    names, d = read_sto(path)
    t = d[:, 0]
    res_f = res_t = 0.0
    saturating = []
    reserves = []
    for j, nm in enumerate(names):
        if j == 0:
            continue
        v = d[:, j]
        if nm in RESIDUAL_FORCE_NAMES:
            res_f = max(res_f, float(np.abs(v).max()))
        elif nm in RESIDUAL_TORQUE_NAMES:
            res_t = max(res_t, float(np.abs(v).max()))
        elif nm in COORD_ACTUATOR_FORCE:
            limit = COORD_ACTUATOR_FORCE[nm] * COORD_ACTUATOR_CONTROL
            peak = float(np.abs(v).max())
            reserves.append((nm, limit, peak, float(v.std())))
            frac_sat = float(np.mean(np.abs(v) >= SATURATION_FRACTION * limit))
            if frac_sat > SATURATION_ALLOWED:
                saturating.append((nm, limit, peak, frac_sat))
    return res_f, res_t, saturating, reserves, t


def check_divergence(path, t0):
    names, d = read_sto(path)
    ang = [j for j, nm in enumerate(names) if j > 0 and "angle" in nm or
           (j > 0 and nm.endswith(("_tilt", "_list", "_rotation", "_flexion_r",
                                   "_adduction_r", "_rotation_r", "_flexion_l",
                                   "_adduction_l", "_rotation_l",
                                   "_extension_r", "_extension_l",
                                   "_dorsiflexion_r", "_dorsiflexion_l",
                                   "_extension", "_bending")))]
    worst, worst_name, has_nan = 0.0, "", False
    for j in ang:
        v = d[:, j]
        if not np.all(np.isfinite(v)):
            has_nan = True
            continue
        span = float(np.abs(v - v[0]).max())
        if span > worst:
            worst, worst_name = span, names[j]
    return worst, worst_name, has_nan


def run_rra_pass(k, n_passes, model_path, t0, t1, ik_window=None, grf_window_xml=None,
                 tasks_xml=None, results_dir=RRA_RESULTS,
                 adjust_com=True):
    """One RRA pass.  Returns (model_for_next_pass, ok, summary)."""
    ik_window = ik_window or IK_WINDOW_MOT
    grf_window_xml = grf_window_xml or GRF_WINDOW_XML
    tasks_xml = tasks_xml or TASKS_XML

    name = "RRA_pass{}".format(k)
    print("\n" + "=" * 78)
    print("RRA pass {}/{}  (model: {})".format(k, n_passes, os.path.basename(model_path)))
    print("=" * 78)

    model = osim.Model(model_path)
    model.updProbeSet().clearAndDestroy()
    model.setName("SQR_RRA")
    model.initSystem()
    if model.getActuators().getSize() == 0:
        # OpenSim writes the adjusted model WITHOUT the residual/reserve force
        # set (verified: 19 actuators in, 0 out), so an adjusted model cannot
        # drive the tracking controller of the next pass.  Fail with the reason
        # instead of letting CMC abort later with "ActuatorForceTarget: ERROR-
        # no controls".  The actuators must be re-applied to it first.
        return None, False, ("{} has no actuators: RRA strips the residual/reserve "
                             "set from its adjusted model, so it cannot seed "
                             "another pass".format(os.path.basename(model_path)))

    tool = osim.RRATool()
    tool.setName(name)
    tool.setModel(model)                       # setModelFilename alone leaves no model
    tool.setDesiredKinematicsFileName(ik_window)
    tool.setTaskSetFileName(tasks_xml)
    # no setConstraintsFileName(): the control limits are on the actuators
    # themselves (min_control/max_control in the model / Actuators.xml)
    tool.setExternalLoadsFileName(grf_window_xml)
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    # -1 = do not filter.  The desired kinematics handed to RRA are already the
    # 15 Hz ik_filter.py output; letting the tool filter them at the same cutoff
    # again would apply the response twice (-6 dB -> -12 dB at 15 Hz), i.e. a
    # different band from the one the GRF is band-limited to.
    tool.setLowpassCutoffFrequency(RRA_LOWPASS_FC)
    tool.setErrorTolerance(ERROR_TOLERANCE)
    tool.setMinDT(MIN_DT)
    tool.setMaxDT(MAX_DT)
    tool.setMaximumNumberOfSteps(int((t1 - t0) * MAX_STEPS_PER_SECOND))
    # RRA only writes an adjusted model when it is allowed to move the COM, and
    # without that file there is nothing to carry into the next pass.  Pass 1
    # therefore always adjusts the COM (as the OpenSim reference does); later
    # passes may disable it once the model has settled.
    tool.setAdjustCOMToReduceResiduals(adjust_com)
    if adjust_com:
        tool.setAdjustedCOMBody(ADJUSTED_COM_BODY)
    tool.setOutputPrecision(8)
    tool.setResultsDir(results_dir)
    # Absolute path: OpenSim resolves a bare file name against the *current
    # working directory*, not against results_dir, so "RRA_pass1_adjusted.osim"
    # used to land next to whoever launched the script and the existence check
    # below never found it -- which silently reduced every run to a single pass.
    adj = os.path.join(results_dir, "{}_adjusted.osim".format(name))
    tool.setOutputModelFileName(adj)

    os.makedirs(results_dir, exist_ok=True)
    try:
        tool.run()
    except Exception as exc:
        return None, False, "RRA threw: {}".format(str(exc).strip().splitlines()[0][:180])

    act = os.path.join(results_dir, "{}_Actuation_force.sto".format(name))
    kin = os.path.join(results_dir, "{}_Kinematics_q.sto".format(name))
    if not os.path.exists(act) or not os.path.exists(kin):
        return None, False, "RRA produced no Actuation_force/Kinematics_q"

    res_f, res_t, sat, reserves, _ = analyse_actuation(act, t0, t1)
    worst, wname, has_nan = check_divergence(kin, t0)

    print("   pelvis residual   : {:.2f} N / {:.2f} Nm   "
          "(target < {:.0f} N / {:.0f} Nm)".format(
              res_f, res_t, RESIDUAL_FORCE_WARN, RESIDUAL_TORQUE_WARN))
    print("   worst coord drift : {:.1f} deg ({}){}".format(
        worst, wname, "  <-- NaN present" if has_nan else ""))
    print("   reserve peaks (N*m/N), cap in brackets; >90% of cap is a problem:")
    for nm, cap, peak, sd in reserves:
        flag = "  <-- NEAR/AT CAP" if peak >= 0.9 * cap else ""
        print("      {:<44s} {:8.2f}  [{:6.1f}]  sd {:7.2f}{}".format(
            nm.split("/")[-1], peak, cap, sd, flag))
    if sat:
        print("   SATURATED actuators:")
        for nm, cap, mx, fr in sat:
            print("      {:<46s} cap {:7.1f}  max|f| {:8.2f}  {:.1f}% of samples".format(
                nm.split("/")[-1], cap, mx, 100 * fr))
    else:
        print("   no reserve actuator saturated")

    if has_nan:
        return None, False, "kinematics contain NaN"
    if worst > DIVERGENCE_DEG:
        return None, False, "diverged: {} moved {:.1f} deg".format(wname, worst)
    if sat:
        return None, False, "{} reserve actuator(s) saturated".format(len(sat))
    if not os.path.exists(adj):
        return None, False, "no adjusted model written (nothing to carry to the next pass)"

    return adj, True, "residual {:.1f} N / {:.1f} Nm, drift {:.1f} deg".format(
        res_f, res_t, worst)


def run_rra_passes(t0, t1, n_passes=DEFAULT_PASSES, ik_window=None, grf_window_xml=None,
                   model_path=None, resume=1, setup_dir=RRA_SETUP_DIR,
                   results_dir=RRA_RESULTS):
    """Run `n_passes` rounds, stopping at the first failure.

    Pass `resume=K` to start at pass K, continuing from
    results_dir/RRA_pass<K-1>_adjusted.osim.  The window files and the setup
    XMLs are looked up in `setup_dir` (what prepare/build_model wrote) unless
    given explicitly.  Returns the history as a list of (pass_number, ok, summary).
    """
    p = paths(setup_dir, results_dir)
    ik_window = ik_window or p["ik_window"]
    grf_window_xml = grf_window_xml or p["grf_window_xml"]

    k_start = max(1, resume)
    if model_path is None:
        if k_start > 1:
            model_path = os.path.join(results_dir,
                                      "RRA_pass{}_adjusted.osim".format(k_start - 1))
            if not os.path.exists(model_path):
                sys.exit("cannot resume: {} not found".format(model_path))
            print("\nResuming at pass {} from {}".format(
                k_start, os.path.basename(model_path)))
        else:
            model_path = p["model"]

    history = []
    for k in range(k_start, n_passes + 1):
        # pass 1 (or the first resumed pass with no adjusted model yet) must be
        # allowed to move the COM, otherwise RRA writes no model to continue from
        adjust_com = (k == k_start)
        nxt, ok, summary = run_rra_pass(k, n_passes, model_path, t0, t1,
                                        ik_window=ik_window,
                                        grf_window_xml=grf_window_xml,
                                        tasks_xml=p["tasks_xml"],
                                        results_dir=results_dir,
                                        adjust_com=adjust_com)
        history.append((k, ok, summary))
        print("   pass {}: {}".format("OK   " if ok else "STOP ", summary))
        if not ok:
            break
        model_path = nxt
    return history


def run(ik_filtered, t0=DEFAULT_TIME_RANGE[0], t1=DEFAULT_TIME_RANGE[1],
        n_passes=DEFAULT_PASSES, resume=1, setup_dir=RRA_SETUP_DIR,
        results_dir=RRA_RESULTS, verbose=True):
    """prepare() + run_rra_passes(): everything RRA needs, from the filtered IK.

    Returns (history, ik_window_path, grf_window_xml_path).
    """
    if verbose:
        print("\nBuilding the RRA model and setup files ...")
    build_model(setup_dir=setup_dir, verbose=verbose)
    grf_win, grf_xml, ik_win = write_window_files(ik_filtered, t0, t1,
                                                  setup_dir=setup_dir)
    if verbose:
        print("\nRunning RRA x{} ...".format(n_passes))
    history = run_rra_passes(t0, t1, n_passes=n_passes, ik_window=ik_win,
                             grf_window_xml=grf_xml, resume=resume,
                             setup_dir=setup_dir, results_dir=results_dir)
    return history, ik_win, grf_xml


def print_summary(history, n_passes=DEFAULT_PASSES, results_dir=RRA_RESULTS):
    print("\n" + "=" * 78)
    print("RRA summary")
    print("=" * 78)
    if n_passes == 0:
        print("\nno RRA passes requested (--passes 0): preparation only")
        print("results in: {}".format(results_dir))
        return
    for k, ok, summary in history:
        print("   pass {} [{}] {}".format(k, "ok" if ok else "FAIL", summary))
    n_ok = sum(1 for _, ok, _ in history if ok)
    print("\n{}/{} passes completed".format(n_ok, len(history)))
    if n_ok == len(history) and history and history[-1][0] == n_passes:
        print("all requested passes finished without divergence or saturation")
    else:
        print("stopped early -- see the last pass above for the reason")
    print("results in: {}".format(results_dir))
