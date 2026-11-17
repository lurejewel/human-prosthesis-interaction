"""Plot: tracked vs no-tracking MocoInverse excitation, 5 sEMG muscles.

Black = processed sEMG reference (e)
Red   = tracked MocoInverse solution (u, 247-parameter tracking run, final)
Blue  = no-tracking MocoInverse baseline (u)
Window: 210-240 s, full 30 s; plus a 3 s zoom panel insets? (single figure, 5 rows)
"""

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
MUSCLES = ["rect_fem_r", "vasti_r", "hamstrings_r", "tib_ant_r", "gastroc_r"]
REF = os.path.join(OUT, "emg_excitation_210_240_right.sto")
TRACKED = os.path.join(OUT, "moco_inverse_solution_B30_cont2_u021_merged.sto")
NOTRACK = os.path.join(OUT, "noemg_baseline", "moco_inverse_solution_nb_u001_merged.sto")
S_TRACKED = [2.0, 1.5084191804226441, 0.9039558238055926, 1.2252940054096981, 1.71818076951833]


def read_sto(path):
    lines = open(path, encoding="utf-8").readlines()
    ih = [i for i, ln in enumerate(lines)
          if ln.strip().lower().startswith("endheader")][0]
    labels = lines[ih + 1].strip().split("\t")
    data = np.genfromtxt(path, skip_header=ih + 2)
    return data[:, 0], {lab: data[:, i] for i, lab in enumerate(labels)}


def main():
    import opensim as osim

    rt, rc = read_sto(REF)
    Tt = osim.MocoTrajectory(TRACKED)
    Tn = osim.MocoTrajectory(NOTRACK)
    tt = np.array(Tt.getTimeMat()).reshape(-1)
    tn = np.array(Tn.getTimeMat()).reshape(-1)

    fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True)
    for ax, m in zip(axes, MUSCLES):
        e = np.interp(tt, rt, rc["/forceset/" + m])
        ut = np.array(Tt.getControlMat("/forceset/" + m)).reshape(-1)
        un = np.array(Tn.getControlMat("/forceset/" + m)).reshape(-1)
        ax.plot(tt, e, color="k", lw=1.6, label="sEMG reference (e)")
        ax.plot(tt, ut, color="tab:red", lw=1.2, alpha=0.85,
                label="tracked MocoInverse (u)")
        ax.plot(tn, un, color="tab:blue", lw=1.2, alpha=0.85,
                label="no-tracking baseline (u)")
        ax.set_ylabel(m, fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)
        if m == MUSCLES[0]:
            ax.legend(loc="upper right", fontsize=8, ncol=3)
    axes[-1].set_xlabel("time (s)")
    axes[-1].set_xlim(210, 240)
    fig.suptitle("B30 30 s right leg: sEMG-driven MocoInverse vs no-tracking baseline", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "..", "..", "results", "moco_tracking_vs_none.png")
    out = os.path.abspath(out)
    fig.savefig(out, dpi=150)
    print("saved", out)


if __name__ == "__main__":
    main()
