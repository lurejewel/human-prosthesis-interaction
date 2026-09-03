"""Quick sanity check for the OpenSim Python API.

Run with the `opensim_scripting` conda environment:
    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe tests\\test_opensim_api.py
"""
import os
import opensim as osim

print("OpenSim version :", osim.GetVersion())
print("Module path     :", osim.__file__)

base_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "experiments", "data", "SQR_walking",
)

for name in ["base_model.osim", "SQR_simbody.osim"]:
    path = os.path.join(base_dir, name)
    print("\n=== Loading:", name, "===")
    model = osim.Model(path)
    model.initSystem()
    print("  Model name     :", model.getName())
    print("  Num bodies     :", model.getNumBodies())
    print("  Num coordinates:", model.getNumCoordinates())
    print("  Num forces     :", model.getForceSet().getSize())
    coords = [
        str(model.getCoordinateSet().get(i).getName())
        for i in range(model.getNumCoordinates())
    ]
    print("  Coordinates    :", ", ".join(coords))
    forces = [
        str(model.getForceSet().get(i).getName())
        for i in range(model.getForceSet().getSize())
    ]
    print("  Forces         :", ", ".join(forces))

print("\nAll OpenSim API checks passed.")
