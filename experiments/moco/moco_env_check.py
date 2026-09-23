"""Environment / data readiness check for OpenSim Moco in this workspace.

Answers one question: *can I run MocoInverse here, and if not, what is missing?*

Run with the ``opensim_scripting`` conda environment::

    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe \\
        experiments\\moco\\moco_env_check.py
    ... --solve        # additionally solve a 1-DOF pendulum with Moco (fast)
    ... --inverse      # additionally run a very short MocoInverse (slow)

Sections
--------
  1. Interpreter + OpenSim + Moco versions
  2. CasADi solver backend (MocoInverse cannot run without it)
  3. Optional Python packages
  4. Input files for the SQR walking dataset
  5. Model <-> inverse-kinematics consistency
  6. Kinematics smoothness *after the example's 15 Hz filter* (the usual cause
     of a failed MocoInverse)
  7. Ground-reaction-force coverage inside the analysis window
  8. ModelProcessor pipeline: actuators added, contact forces present
  9. A real (tiny) optimal-control solve, to prove the solver stack works

Sections 6, 8 and 9b are supposed to predict the *real* run, so they do not
re-declare the example script's numbers: the residual capacities, the
kinematics filter cutoff, the analysis window and the radians conversion are
all taken from ``moco_inverse.py`` (imported below).  Using the OpenSim
tutorial values here instead is what previously made this check disagree with
the actual solve.

Exit code is 0 when every CRITICAL check passes, 1 otherwise.
"""

import argparse
import os
import sys
import time
import traceback

import numpy as np

# Locate the CasADi solver plugins before importing opensim, so that the check
# works even when the conda environment was not activated.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moco_bootstrap                                    # noqa: E402
CASADI_PLUGIN_DIR = moco_bootstrap.ensure_casadi_plugins()

import opensim as osim                                   # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "experiments", "data", "SQR_walking")

MODEL_FILE = os.path.join(DATA_DIR, "SQR_simbody.osim")
# 15 Hz filtered kinematics (ik_filter.py); level_walking_ik.mot is the raw one.
IK_FILE = os.path.join(DATA_DIR, "level_walking_ik_filtered.mot")
GRF_MOT = os.path.join(DATA_DIR, "level_walking_grf_filtered.mot")
GRF_XML = os.path.join(DATA_DIR, "level_walking_grf_filtered.xml")

# ---------------------------------------------------------------------------
# Parameters shared with the script under test
# ---------------------------------------------------------------------------
# This file exists to predict whether ``moco_inverse.py`` will run, so
# the model it builds and the preprocessing it applies must match that script
# exactly.  Re-declaring the numbers here is what previously made the two
# disagree:
#   * the residual capacities were the OpenSim tutorial values (250 N*m / 50 N)
#     while the example deliberately uses larger ones -- with 50 N the pelvis
#     translational residuals saturate and the problem is infeasible, so a
#     failure reported here said nothing about the example script;
#   * the smoke MocoInverse was fed the raw, degree-based IK file, without the
#     radians conversion and without the 15 Hz low-pass the example applies,
#     either of which on its own prevents convergence.
# Importing the parameters (and reusing the example's helper functions) keeps
# the two files from drifting apart again.
try:
    import moco_inverse as example                       # noqa: E402
    RESIDUAL_ROT         = example.RESIDUAL_ROTATIONAL_FORCE
    RESIDUAL_TRANS       = example.RESIDUAL_TRANSLATIONAL_FORCE
    KINEMATICS_FILTER_HZ = example.KINEMATICS_FILTER_HZ
    _example_defaults    = example.parse_args([])                # defaults only
    WINDOW               = (_example_defaults.start, _example_defaults.end)
    EXAMPLE_IMPORT_ERROR = None
except Exception as _exc:                   # scipy missing, syntax error, ...
    example              = None
    RESIDUAL_ROT         = 400.0            # keep in sync with the example
    RESIDUAL_TRANS       = 500.0            # keep in sync with the example
    KINEMATICS_FILTER_HZ = 15.0             # keep in sync with the example
    WINDOW               = (180.0, 181.5)   # keep in sync with the example
    EXAMPLE_IMPORT_ERROR = _exc

CRITICAL, WARN, INFO = "CRITICAL", "WARN", "INFO"
_findings = []


def record(level, section, message):
    _findings.append((level, section, message))
    marker = {CRITICAL: "[FAIL]", WARN: "[warn]", INFO: "[ ok ]"}[level]
    print(f"{marker} {section:<34} {message}")


def section(title):
    print(f"\n{'-' * 78}\n{title}\n{'-' * 78}")


# ---------------------------------------------------------------------------
def check_versions():
    section("1. Interpreter and libraries")
    record(INFO, "Python", f"{sys.version.split()[0]}  ({sys.executable})")
    record(INFO, "OpenSim", osim.GetVersion())
    record(INFO, "Moco", osim.GetMocoVersion())
    if sys.version_info[:2] >= (3, 13):
        record(WARN, "Python version",
               "conda OpenSim packages target 3.9-3.12; 3.13+ may not be "
               "installable")


def check_casadi():
    section("2. CasADi solver backend")
    try:
        available = osim.MocoCasADiSolver.isAvailable()
    except Exception as exc:
        record(CRITICAL, "MocoCasADiSolver", f"probe raised: {exc}")
        return
    if available:
        record(INFO, "MocoCasADiSolver.isAvailable()", "True")
    else:
        record(CRITICAL, "MocoCasADiSolver.isAvailable()",
               "False -- MocoInverse/MocoTrack CANNOT run. Install the "
               "opensim package from the opensim-org conda channel.")
        return

    # isAvailable() is a compile-time flag: it stays True even when the CasADi
    # *plugin* libraries (IPOPT) cannot be loaded at runtime, which is exactly
    # what happens when python.exe is run without `conda activate`.
    if CASADI_PLUGIN_DIR:
        record(INFO, "CasADi plugin dir", CASADI_PLUGIN_DIR)
    else:
        record(CRITICAL, "CasADi plugin dir",
               "casadi_nlpsol_ipopt.dll not found -- every solve will fail "
               "with \"Plugin 'ipopt' is not found\". Run "
               "'conda activate opensim_scripting' or set CASADIPATH.")

    try:
        solver = osim.MocoCasADiSolver()
        record(INFO, "instantiate MocoCasADiSolver", solver.getClassName())
    except Exception as exc:
        record(CRITICAL, "instantiate MocoCasADiSolver", str(exc))

    # The CasADi *runtime* is bundled as DLLs; the python 'casadi' module is a
    # separate thing and is NOT required by Moco.
    try:
        import casadi  # noqa: F401
        record(INFO, "python 'casadi' module", f"{casadi.__version__} (optional)")
    except ImportError:
        record(INFO, "python 'casadi' module",
               "not installed -- not required by Moco (only needed by "
               "plot_casadi_sparsity.py)")


def check_packages():
    section("3. Optional Python packages")
    for mod, why in [
        ("numpy", "required by the example script"),
        ("scipy", "used by run_ik_id.py filtering"),
        ("matplotlib", "required by osim.report / plotting"),
        ("pandas", "convenient for result tables"),
        ("sklearn", "only for the muscle-synergy example"),
    ]:
        try:
            m = __import__(mod)
            record(INFO, mod, f"{getattr(m, '__version__', '?')}  ({why})")
        except ImportError:
            lvl = CRITICAL if mod in ("numpy",) else WARN
            record(lvl, mod, f"MISSING ({why})")


def check_files():
    section("4. Input files (SQR walking dataset)")
    for label, path in [("model", MODEL_FILE), ("IK kinematics", IK_FILE),
                        ("GRF (15 Hz)", GRF_MOT), ("external loads xml", GRF_XML)]:
        if os.path.exists(path):
            record(INFO, label, f"{path}  ({os.path.getsize(path) / 1024:.0f} KB)")
        else:
            record(CRITICAL, label, f"MISSING: {path}")
    return all(os.path.exists(p) for p in (MODEL_FILE, IK_FILE, GRF_XML))


def check_model_and_ik():
    section("5. Model <-> kinematics consistency")
    model = osim.Model(MODEL_FILE)
    model.initSystem()

    coords = [model.getCoordinateSet().get(i).getName()
              for i in range(model.getCoordinateSet().getSize())]
    muscles = [model.getMuscles().get(i).getName()
               for i in range(model.getMuscles().getSize())]

    table = osim.TimeSeriesTable(IK_FILE)
    ik_cols = list(table.getColumnLabels())
    times = np.array(table.getIndependentColumn())

    record(INFO, "model coordinates", str(len(coords)))
    record(INFO, "model muscles", str(len(muscles)))
    record(INFO, "IK columns", str(len(ik_cols)))
    record(INFO, "IK time span", f"{times[0]:.3f} - {times[-1]:.3f} s, "
                                 f"{len(times)} rows "
                                 f"({1.0 / np.median(np.diff(times)):.0f} Hz)")

    missing = sorted(set(coords) - set(ik_cols))
    extra = sorted(set(ik_cols) - set(coords))
    if missing:
        record(CRITICAL, "coordinates absent from IK file",
               f"{missing} -- MocoInverse requires ALL coordinates")
    else:
        record(INFO, "all coordinates present in IK file", "yes")
    if extra:
        record(WARN, "extra IK columns",
               f"{extra} -- fine, set_kinematics_allow_extra_columns(True)")
    else:
        record(INFO, "extra IK columns", "none")

    # NaN check -- MocoInverse cannot handle NaNs in the kinematics.
    data = table.getMatrix().to_numpy()
    if np.isnan(data).any():
        record(CRITICAL, "NaNs in IK data", f"{int(np.isnan(data).sum())} values")
    else:
        record(INFO, "NaNs in IK data", "none")

    return model, coords, ik_cols, table


def lowpass_columns(q, t, cutoff_hz):
    """Zero-phase 4th-order Butterworth low-pass, exactly as the example does.

    ``q`` is the *whole* table (every row), because the example filters the
    entire IK record and only then restricts to the analysis window; filtering
    a window-sized slice instead would add edge transients the solver never
    sees.  Returns ``(filtered, note)`` where ``note`` is meant for the log.
    """
    if np.isnan(q).any():
        return q, "NaNs in the IK data -- filtering skipped (see section 5)"
    if not cutoff_hz or cutoff_hz <= 0:
        return q, "filtering disabled (example --filter-hz 0)"
    fs = 1.0 / float(np.median(np.diff(t)))
    if cutoff_hz >= fs / 2.0:
        return q, (f"cutoff {cutoff_hz} Hz >= Nyquist {fs / 2:.0f} Hz -- the "
                   f"example skips the filter too")
    try:
        from scipy.signal import butter, filtfilt
    except ImportError:
        return q, "scipy NOT installed -- unfiltered, so the values below are raw"
    b, a = butter(4, cutoff_hz / (fs / 2.0), btype="low")
    return filtfilt(b, a, q, axis=0), f"{cutoff_hz} Hz zero-phase (fs={fs:.0f} Hz)"


def check_smoothness(table, coords):
    section("6. Kinematics smoothness (accelerations MocoInverse must reproduce)")
    times = np.array(table.getIndependentColumn())
    data = table.getMatrix().to_numpy()
    labels = list(table.getColumnLabels())

    # MocoInverse differentiates the kinematics twice, and the example script
    # therefore low-passes the coordinates at KINEMATICS_FILTER_HZ before
    # handing them over.  Report the accelerations *after* that filter, i.e.
    # what the solver actually sees, and keep the raw value alongside it so the
    # effect of the filter is visible.
    data_f, note = lowpass_columns(data, times, KINEMATICS_FILTER_HZ)
    record(INFO, "preprocessing (as in the example)", note)

    # Restrict to the analysis window.
    sel = (times >= WINDOW[0]) & (times <= WINDOW[1])
    if sel.sum() < 5:
        record(WARN, "analysis window", f"only {int(sel.sum())} samples in "
                                        f"{WINDOW[0]}-{WINDOW[1]} s")
        return
    t = times[sel]
    dt = float(np.median(np.diff(t)))

    acc = np.gradient(np.gradient(data_f[sel, :], dt, axis=0), dt, axis=0)
    peak = np.max(np.abs(acc), axis=0)
    acc_raw = np.gradient(np.gradient(data[sel, :], dt, axis=0), dt, axis=0)
    peak_raw = np.max(np.abs(acc_raw), axis=0)

    order = np.argsort(-peak)
    record(INFO, "peak |q''| in window [deg/s^2 or m/s^2]",
           f"max={peak[order[0]]:.1f} on '{labels[order[0]]}' "
           f"(raw, unfiltered: {peak_raw[order[0]]:.1f})")
    hot = [(labels[i], peak[i], peak_raw[i]) for i in order[:5]]
    for name, val, val_raw in hot:
        lvl = WARN if val > 5000 else INFO
        record(lvl, f"  {name}", f"peak |q''| = {val:.1f} (raw {val_raw:.1f})")

    if peak[order[0]] > 20000:
        record(WARN, "smoothness verdict",
               f"very large accelerations even after the {KINEMATICS_FILTER_HZ} "
               f"Hz filter -- lower --filter-hz or improve the IK "
               f"(MocoInverse differentiates the kinematics twice)")
    else:
        record(INFO, "smoothness verdict",
               f"acceptable after the {KINEMATICS_FILTER_HZ} Hz filter "
               f"(this is what MocoInverse will differentiate)")


def check_grf():
    section("7. Ground reaction force coverage in the window")
    if not os.path.exists(GRF_MOT):
        record(WARN, "GRF file", "missing, skipped")
        return
    tbl = osim.TimeSeriesTable(GRF_MOT)
    t = np.array(tbl.getIndependentColumn())
    labels = list(tbl.getColumnLabels())
    sel = (t >= WINDOW[0]) & (t <= WINDOW[1])
    if not sel.any():
        record(WARN, "window", f"no GRF samples in {WINDOW}")
        return

    for tag in ("ground_force1_vy", "ground_force2_vy"):
        if tag not in labels:
            record(WARN, tag, "column not found")
            continue
        col = tbl.getDependentColumn(tag).to_numpy()[sel]
        contact = float(np.mean(col > 10.0)) * 100.0
        record(INFO, tag, f"peak={col.max():.1f} N, "
                          f"in-contact {contact:.0f}% of window")
    record(INFO, "window checked", f"{WINDOW[0]}-{WINDOW[1]} s")


def build_example_processor(grf_xml):
    """The ModelProcessor pipeline of ``moco_inverse.py``.

    Reused by section 8 and by the ``--inverse`` smoke solve, so that this
    check, the smoke solve and the real run all build the same model.  When the
    example module could not be imported, the same operator chain is rebuilt
    here with the shared residual capacities (the values that matter).
    """
    if example is not None:
        return example.build_model_processor(grf_xml,
                                             residual_rot=RESIDUAL_ROT,
                                             residual_trans=RESIDUAL_TRANS)

    mp = osim.ModelProcessor(MODEL_FILE)
    mp.append(osim.ModOpAddExternalLoads(grf_xml))
    mp.append(osim.ModOpIgnoreTendonCompliance())
    mp.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    mp.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    mp.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    mp.append(osim.ModOpAddResiduals(RESIDUAL_ROT, RESIDUAL_TRANS, 1.0))
    mp.append(osim.ModOpAddReserves(1.0))
    return mp


def check_model_processing():
    section("8. ModelProcessor pipeline (as used by the example script)")
    mp = build_example_processor(GRF_XML)
    record(INFO, "pipeline source",
           "moco_inverse.build_model_processor()" if example is not None
           else f"local rebuild (import failed: {EXAMPLE_IMPORT_ERROR})")
    record(INFO, "residual capacity",
           f"{RESIDUAL_ROT:g} N*m / {RESIDUAL_TRANS:g} N (pelvis)")

    try:
        model = mp.process()
        model.initSystem()
    except Exception as exc:
        record(CRITICAL, "ModelProcessor.process()", str(exc))
        return None

    n_musc = model.getMuscles().getSize()
    names = [model.getActuators().get(i).getName()
             for i in range(model.getActuators().getSize())]
    residuals = [n for n in names if n.startswith("residual_")]
    reserves = [n for n in names if n.startswith("reserve_")]

    record(INFO, "muscles after conversion", f"{n_musc} "
           f"({model.getMuscles().get(0).getConcreteClassName()})")
    record(INFO, "residual actuators", f"{len(residuals)} (pelvis DOF)")
    record(INFO, "reserve actuators", f"{len(reserves)} (other DOF)")
    record(INFO, "total actuators", str(len(names)))

    ncoord = model.getNumCoordinates()
    if len(residuals) + len(reserves) >= ncoord:
        record(INFO, "every coordinate actuated", f"{ncoord} coords covered")
    else:
        record(WARN, "coordinate coverage",
               f"{len(residuals) + len(reserves)} coordinate actuators for "
               f"{ncoord} coordinates -- some DOF rely on muscles only")

    contacts = [model.getForceSet().get(i).getName()
                for i in range(model.getForceSet().getSize())
                if "HuntCrossley" in model.getForceSet().get(i).getConcreteClassName()]
    if contacts:
        record(WARN, "contact forces in model",
               f"{contacts} -- DISABLE them when prescribing measured GRF, "
               f"otherwise the ground force is counted twice")
    else:
        record(INFO, "contact forces in model", "none")

    limits = [model.getForceSet().get(i).getName()
              for i in range(model.getForceSet().getSize())
              if "CoordinateLimitForce" in model.getForceSet().get(i).getConcreteClassName()]
    record(INFO, "coordinate limit forces", str(limits or "none"))
    return model


def solve_pendulum():
    """Solve a trivial 1-DOF optimal control problem end to end.

    This exercises CasADi + IPOPT without touching the user's data, in ~1 s.
    """
    section("9a. Solver smoke test (1-DOF pendulum, ~1 s)")
    try:
        model = osim.ModelFactory.createPendulum()
        model.setName("pendulum")

        study = osim.MocoStudy()
        study.setName("pendulum_smoke")
        problem = study.updProblem()
        problem.setModel(model)
        problem.setTimeBounds(0.0, 1.0)
        problem.setStateInfo("/jointset/j0/q0/value", [-10, 10], 0.0, 0.0)
        problem.setStateInfo("/jointset/j0/q0/speed", [-50, 50], 0.0, 0.0)
        problem.setControlInfo("/tau0", [-100, 100])
        problem.addGoal(osim.MocoControlGoal("control_effort", 1.0))

        # updSolver() is typed as the MocoSolver base class in Python, so the
        # CasADi-specific settings need an explicit downcast.
        study.initCasADiSolver()
        solver = osim.MocoCasADiSolver.safeDownCast(study.updSolver())
        solver.set_optim_max_iterations(200)

        t0 = time.time()
        solution = study.solve()
        elapsed = time.time() - t0
        if solution.success():
            record(INFO, "pendulum solve", f"SUCCESS in {elapsed:.1f} s, "
                                           f"objective={solution.getObjective():.4g}")
            return True
        record(CRITICAL, "pendulum solve", "solver returned unsuccessful")
        return False
    except Exception as exc:
        record(CRITICAL, "pendulum solve", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
        return False


def run_tiny_inverse():
    """Run a very short MocoInverse on the real data (slow, optional).

    This mirrors ``moco_inverse.py`` rather than approximating it: the
    same ModelProcessor, the same contact forces disabled, and the same
    kinematics preprocessing (whole-file 15 Hz low-pass + deg -> rad).  Feeding
    the raw degree-based IK file straight to ``setKinematics`` -- as this
    function used to -- asks the solver for a 67 radian knee angle, which
    cannot converge, so a failure here said nothing about the real script.
    """
    section("9b. Real MocoInverse smoke test (0.3 s window)")
    if example is None:
        record(WARN, "MocoInverse smoke test",
               f"skipped: moco_inverse.py could not be imported "
               f"({EXAMPLE_IMPORT_ERROR}); this check needs its preprocessing")
        return False
    try:
        model = build_example_processor(GRF_XML).process()
        disabled = example.disable_contact_forces(model)
        model.initSystem()
        record(INFO, "contacts disabled", str(disabled or "(none found)"))

        # The example's own preprocessing, written one directory over so it
        # cannot clobber the example's output folder.  (The helper does not
        # create its target directory -- the example's main() does that.)
        out_dir = os.path.join(SCRIPT_DIR, "outputs", "env_check")
        os.makedirs(out_dir, exist_ok=True)
        ik_file = example.prepare_kinematics_radians(
            model, IK_FILE, out_dir, filter_hz=KINEMATICS_FILTER_HZ,
            tag="ik_radians")

        inverse = osim.MocoInverse()
        inverse.setName("tiny_inverse")
        inverse.setModel(osim.ModelProcessor(model))
        inverse.setKinematics(osim.TableProcessor(ik_file))
        inverse.set_initial_time(WINDOW[0])
        inverse.set_final_time(WINDOW[0] + 0.3)
        inverse.set_mesh_interval(0.05)
        inverse.set_kinematics_allow_extra_columns(True)
        inverse.set_max_iterations(500)

        t0 = time.time()
        solution = inverse.solve()
        elapsed = time.time() - t0
        sol = solution.getMocoSolution()
        if sol.success():
            record(INFO, "MocoInverse solve", f"SUCCESS in {elapsed:.1f} s, "
                                              f"objective={sol.getObjective():.4g}")
            return True
        record(WARN, "MocoInverse solve",
               f"finished without reaching tolerance in {elapsed:.1f} s. This "
               f"runs the example's own pipeline (0.3 s window, 50 ms mesh), so "
               f"the failure is meaningful -- see the solver log and compare "
               f"with the full example run before drawing conclusions")
        return False
    except Exception as exc:
        record(CRITICAL, "MocoInverse solve", f"{type(exc).__name__}: {exc}")
        traceback.print_exc()
        return False


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description="Check whether Moco can run here.")
    p.add_argument("--solve", action="store_true",
                   help="also solve a trivial pendulum problem")
    p.add_argument("--inverse", action="store_true",
                   help="also run a very short MocoInverse on the real data")
    args = p.parse_args(argv)

    print("=" * 78)
    print("OpenSim Moco environment check")
    print("=" * 78)

    check_versions()
    check_casadi()
    check_packages()

    files_ok = check_files()
    model = coords = table = None
    if files_ok:
        try:
            model, coords, ik_cols, table = check_model_and_ik()
        except Exception as exc:
            record(CRITICAL, "model/IK load", f"{type(exc).__name__}: {exc}")
        if table is not None:
            check_smoothness(table, coords)
        check_grf()
    check_model_processing()

    solved = None
    if args.solve or args.inverse:
        solved = solve_pendulum()
    if args.inverse:
        run_tiny_inverse()

    # ---- verdict ----------------------------------------------------------
    section("VERDICT")
    criticals = [f for f in _findings if f[0] == CRITICAL]
    warns = [f for f in _findings if f[0] == WARN]

    if criticals:
        print("Moco CANNOT run as configured. Blocking problems:")
        for _, sec, msg in criticals:
            print(f"  * {sec}: {msg}")
        print("\nSee experiments/moco/README.md for installation instructions.")
        return 1

    print("Moco IS available and the inputs look usable.")
    print(f"  {len(warns)} warning(s):")
    for _, sec, msg in warns:
        print(f"  * {sec}: {msg}")
    if solved is False:
        print("  * the optional solve did not succeed (see above)")
    print("\nNext step:")
    print("  python experiments/moco/moco_inverse.py --smoke")
    return 0


if __name__ == "__main__":
    sys.exit(main())
