"""Make OpenSim Moco usable without an activated conda environment.

The problem
-----------
The conda ``opensim`` package ships CasADi as shared libraries under
``<env>/Library/bin`` -- including the nonlinear-solver plugins such as
``casadi_nlpsol_ipopt.dll``.  ``conda activate`` puts that directory on PATH,
so everything works from an activated shell.  But if ``python.exe`` is launched
directly (VS Code "Run Python File", a task runner, ``subprocess``, a batch
job, ...) that directory is *not* on PATH and the failure is confusing:

* ``opensim.MocoCasADiSolver.isAvailable()`` still returns ``True``, because
  that flag is set at compile time, and
* ``osim.MocoCasADiSolver()`` still constructs fine, but
* every actual solve dies with::

      MocoCasADiSolver failed internally with message:
      Plugin 'ipopt' is not found.

``ensure_casadi_plugins()`` locates the plugin directory next to the running
interpreter and exports it through ``CASADIPATH``, ``PATH`` and
``os.add_dll_directory()``.  Call it *before* ``import opensim``.

Typical use::

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import moco_bootstrap
    moco_bootstrap.ensure_casadi_plugins()
    import opensim as osim
"""

import os
import sys

PLUGIN_NAME = "casadi_nlpsol_ipopt.dll"
_dll_handles = []          # keep add_dll_directory handles alive for the process


def _candidate_dirs():
    """Directories that could hold the CasADi plugin libraries, best first."""
    exe_dir = os.path.dirname(os.path.abspath(sys.executable))
    env_dirs = []
    for root in (exe_dir, sys.prefix, os.path.dirname(exe_dir)):
        if root and root not in env_dirs:
            env_dirs.append(root)

    dirs = []
    if os.environ.get("CASADIPATH"):
        dirs.append(os.environ["CASADIPATH"])
    for root in env_dirs:
        dirs.append(os.path.join(root, "Library", "bin"))
        dirs.append(os.path.join(root, "Library", "casadi", "casadi"))
        dirs.append(os.path.join(root, "bin"))
        dirs.append(os.path.join(root, "lib"))
    # Anything already on PATH (the `conda activate` case).
    dirs.extend(p for p in os.environ.get("PATH", "").split(os.pathsep) if p)
    return dirs


def find_casadi_plugin_dir():
    """Return the first directory containing the CasADi IPOPT plugin, else None."""
    seen = set()
    for d in _candidate_dirs():
        key = os.path.normcase(os.path.abspath(d))
        if key in seen:
            continue
        seen.add(key)
        try:
            if os.path.isfile(os.path.join(d, PLUGIN_NAME)):
                return os.path.abspath(d)
        except OSError:
            continue
    return None


def ensure_casadi_plugins(verbose=False):
    """Point CasADi at its plugin directory.

    Returns the directory that was registered, or ``None`` when it could not be
    found (in which case the caller should treat solves as unavailable).
    """
    plugin_dir = find_casadi_plugin_dir()
    if plugin_dir is None:
        if verbose:
            print("[moco_bootstrap] CasADi plugin directory NOT found "
                  f"(looked for {PLUGIN_NAME}); solves will fail with "
                  "'Plugin ipopt is not found'. Activate the conda env or set "
                  "CASADIPATH.", file=sys.stderr)
        return None

    # CASADIPATH is the documented search location; PATH covers the loader's
    # fallback search; add_dll_directory helps on Python 3.8+ Windows.
    if os.environ.get("CASADIPATH") != plugin_dir:
        os.environ["CASADIPATH"] = plugin_dir

    path_parts = os.environ.get("PATH", "").split(os.pathsep)
    if not any(os.path.normcase(p) == os.path.normcase(plugin_dir)
               for p in path_parts if p):
        os.environ["PATH"] = plugin_dir + os.pathsep + os.environ.get("PATH", "")

    if hasattr(os, "add_dll_directory"):
        try:
            _dll_handles.append(os.add_dll_directory(plugin_dir))
        except (OSError, AttributeError):
            pass

    if verbose:
        print(f"[moco_bootstrap] CasADi plugins: {plugin_dir}")
    return plugin_dir
