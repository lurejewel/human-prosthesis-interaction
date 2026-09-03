"""Temporary GRF spectral analysis for vx..tz (belt1 + belt2), interactive plots.

Shows an interactive matplotlib figure window (zoom / pan enabled) — no PNG is
saved by default.  Run with the `opensim_scripting` conda env:

    D:\\Software\\miniconda_py312\\envs\\opensim_scripting\\python.exe "experiments\\grf_spectrum_analysis.py"

The .mot file has 18 GRF data columns (see run_ik_id.py for the layout):
  belt1: ground_force1_{vx,vy,vz,px,py,pz}, ground_torque1_{x,y,z}
  belt2: ground_force2_{vx,vy,vz,px,py,pz}, ground_torque2_{x,y,z}
"""
import argparse
import os

import numpy as np
from scipy.signal import detrend, find_peaks, welch

# ---------------------------------------------------------------------------
# Pick an interactive backend BEFORE importing pyplot so the figure window
# has the zoom / pan toolbar.  TkAgg is usually present in the conda env.
# ---------------------------------------------------------------------------
import matplotlib

_BACKEND = None
for _b in ("TkAgg", "QtAgg", "Qt5Agg"):
    try:
        matplotlib.use(_b)
        _BACKEND = _b
        break
    except Exception:
        continue

import matplotlib.pyplot as plt  # noqa: E402

if _BACKEND:
    print("[backend] interactive matplotlib backend :", _BACKEND)
else:
    print("[warn] no interactive backend found; the window may not have a zoom toolbar")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, "data", "SQR_walking")
GRF_MOT    = os.path.join(DATA_DIR, "level_walking_grf.mot")

COMPONENTS = ["vx (N)", "vy (N)", "vz (N)",
              "px (m)", "py (m)", "pz (m)",
              "tx (Nm)", "ty (Nm)", "tz (Nm)"]
BELT_NAMES = ["belt1", "belt2"]

DEFAULT_T0 = 100.0   # steady-state analysis window (walking starts ~6.3 s)
DEFAULT_T1 = 240.0
DEFAULT_FMAX = 15.0  # default x-axis limit in Hz (you can zoom further)


# ---------------------------------------------------------------------------
# Helpers (read_table / column_names copied from run_ik_id.py)
# ---------------------------------------------------------------------------
def read_table(path):
    """Read an OpenSim .mot/.sto file -> (header_lines, column_labels, data)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    header = []
    i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("endheader"):
        header.append(lines[i].rstrip("\r\n"))
        i += 1
    if i < len(lines):
        header.append("endheader")
        i += 1

    col_labels = None
    while i < len(lines):
        ln = lines[i]
        if ln.strip():
            try:
                float(ln.lstrip().split()[0])
            except ValueError:
                col_labels = ln.strip()
                i += 1
            break
        i += 1

    data_lines = [ln for ln in lines[i:] if ln.strip()]
    data = np.loadtxt(data_lines)
    return header, col_labels, data


def column_names(col_labels):
    if not col_labels:
        return []
    return [c.strip() for c in col_labels.split()]


def get_component_columns(cols):
    """Map column names -> (belt1_idx, belt2_idx), each 9 ints in
    [vx,vy,vz,px,py,pz,tx,ty,tz] order.  Returns None if names are missing."""
    if not cols or len(cols) < 19:
        return None
    idx = {name: i for i, name in enumerate(cols)}
    needed1 = (["ground_force1_" + s for s in ("vx", "vy", "vz", "px", "py", "pz")]
               + ["ground_torque1_" + s for s in ("x", "y", "z")])
    needed2 = (["ground_force2_" + s for s in ("vx", "vy", "vz", "px", "py", "pz")]
               + ["ground_torque2_" + s for s in ("x", "y", "z")])
    if all(n in idx for n in needed1 + needed2):
        return ([idx[n] for n in needed1], [idx[n] for n in needed2])
    return None


# ---------------------------------------------------------------------------
# Spectral analysis
# ---------------------------------------------------------------------------
def single_sided_spectrum(x, fs):
    """Amplitude spectrum (single-sided), DC removed + linear detrend."""
    x = x - np.mean(x)
    x = detrend(x)
    n = len(x)
    xf = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    amp = np.abs(xf) / n
    amp[1:] *= 2.0            # single-sided (double everything but DC)
    return freqs, amp


def detect_peaks(freqs, amp, fmax, min_freq=0.2, topn=5):
    """Return the top `topn` peaks in [min_freq, fmax] as [(f, a), ...]."""
    sel = (freqs >= min_freq) & (freqs <= fmax)
    fr, a = freqs[sel], amp[sel]
    if len(a) == 0:
        return []
    df = freqs[1] - freqs[0]
    height = 0.05 * a.max()                       # 5% of local max
    distance = max(int(0.4 / df), 2)              # peaks >= 0.4 Hz apart
    pk, _ = find_peaks(a, height=height, distance=distance)
    if len(pk) == 0:
        return []
    order = np.argsort(a[pk])[::-1][:topn]
    return [(float(fr[i]), float(a[i])) for i in pk[order]]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def make_fft_figure(freqs, spectra, peaks_by_belt, fmax, ylog, t0, t1, fs, path):
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    colors = ["C0", "C3"]
    for k, comp in enumerate(COMPONENTS):
        ax = axes[k // 3][k % 3]
        for b in range(2):
            ax.plot(freqs, spectra[b][k], color=colors[b], lw=0.9, label=BELT_NAMES[b])
            for f0, a0 in peaks_by_belt[b][k]:
                ax.annotate("{:.2f} Hz".format(f0), xy=(f0, a0),
                            xytext=(f0, a0 * 1.18), ha="center",
                            fontsize=7, color=colors[b])
        ax.set_title(comp, fontsize=10)
        ax.set_xlim(0, fmax)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
        if ylog:
            ax.set_yscale("log")
    for ax in axes[-1]:
        ax.set_xlabel("frequency (Hz)")
    for ax in axes[:, 0]:
        ax.set_ylabel("amplitude")
    fig.suptitle("GRF amplitude spectrum — {} ({}–{} s, fs={:.0f} Hz)".format(
        os.path.basename(path), t0, t1, fs))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def make_welch_figure(psd_by_belt, fmax, t0, t1, fs, path):
    fig, axes = plt.subplots(3, 3, figsize=(15, 10), sharex=True)
    colors = ["C0", "C3"]
    for k, comp in enumerate(COMPONENTS):
        ax = axes[k // 3][k % 3]
        for b in range(2):
            f, p = psd_by_belt[b][k]
            ax.plot(f, 10.0 * np.log10(p + 1e-30), color=colors[b], lw=0.9,
                    label=BELT_NAMES[b])
        ax.set_title(comp, fontsize=10)
        ax.set_xlim(0, fmax)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper right")
    for ax in axes[-1]:
        ax.set_xlabel("frequency (Hz)")
    for ax in axes[:, 0]:
        ax.set_ylabel("PSD (dB/Hz)")
    fig.suptitle("GRF Welch PSD — {} ({}–{} s, fs={:.0f} Hz)".format(
        os.path.basename(path), t0, t1, fs))
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="GRF spectral analysis (interactive, no PNG)")
    parser.add_argument("--file", default=GRF_MOT, help="GRF .mot file (default: raw level_walking_grf.mot)")
    parser.add_argument("--belt", choices=["1", "2", "both"], default="both",
                        help="which belt(s) to plot")
    parser.add_argument("--t0", type=float, default=DEFAULT_T0, help="window start [s]")
    parser.add_argument("--t1", type=float, default=DEFAULT_T1, help="window end [s]")
    parser.add_argument("--full", action="store_true", help="use the whole time series")
    parser.add_argument("--fmax", type=float, default=DEFAULT_FMAX, help="x-axis limit in Hz")
    parser.add_argument("--ylog", action="store_true", help="log y-axis")
    parser.add_argument("--welch", action="store_true", help="also open a Welch PSD figure")
    parser.add_argument("--save", default=None, help="optional fallback: save figure to PNG")
    parser.add_argument("--no-show", action="store_true", help="compute + print only (no window)")
    args = parser.parse_args()

    # --- read data ---------------------------------------------------------
    header, col_labels, data = read_table(args.file)
    t = data[:, 0]
    fs = 1.0 / np.median(np.diff(t))

    cols = column_names(col_labels)
    mapping = get_component_columns(cols)
    if mapping is not None:
        c1, c2 = mapping
    else:                                    # positional fallback (belt1=1..9, belt2=10..18)
        c1 = list(range(1, 10))
        c2 = list(range(10, 19))
        print("[warn] column names not recognized, using positional mapping")

    # --- analysis window ---------------------------------------------------
    if args.full:
        t0, t1 = t[0], t[-1]
    else:
        t0, t1 = args.t0, args.t1
    mask = (t >= t0) & (t <= t1)
    if mask.sum() < 100:
        raise ValueError("analysis window [{} {}] contains too few samples".format(t0, t1))
    print("data   : {} rows, fs = {:.0f} Hz, window [{:.2f} {:.2f}] s ({} samples)"
          .format(len(t), fs, t0, t1, mask.sum()))
    print("resolve: df = {:.4f} Hz  (1/window_duration)".format(1.0 / (t1 - t0)))

    belts = {"1": [0], "2": [1], "both": [0, 1]}[args.belt]
    comp_idx = [c1, c2]

    # --- spectra -----------------------------------------------------------
    spectra = [[None] * 9 for _ in range(2)]     # [belt][comp]
    psds    = [[None] * 9 for _ in range(2)]
    peaks   = [[None] * 9 for _ in range(2)]
    freqs = None

    for b in belts:
        for k in range(9):
            x = data[mask, comp_idx[b][k]]
            freqs, amp = single_sided_spectrum(x, fs)
            spectra[b][k] = amp
            peaks[b][k] = detect_peaks(freqs, amp, args.fmax)
            if args.welch:
                f, p = welch(x, fs=fs, nperseg=min(2 ** 12, len(x) // 4),
                             detrend="linear", scaling="density")
                psds[b][k] = (f, p)

    # --- print dominant frequencies ----------------------------------------
    print("\nDominant peaks (top 5, {}-{} Hz):".format(0.2, args.fmax))
    for k, comp in enumerate(COMPONENTS):
        parts = []
        for b in belts:
            desc = ", ".join("{:.2f} Hz({:.3g})".format(f0, a0)
                             for f0, a0 in peaks[b][k]) or "—"
            parts.append("{}: {}".format(BELT_NAMES[b], desc))
        print("  {:<9} | ".format(comp) + " | ".join(parts))

    # --- figures ------------------------------------------------------------
    figs = []
    figs.append(make_fft_figure(freqs, spectra, peaks, args.fmax, args.ylog,
                                t0, t1, fs, args.file))
    if args.welch:
        figs.append(make_welch_figure(psds, args.fmax, t0, t1, fs, args.file))

    if args.save:
        os.makedirs(os.path.dirname(os.path.abspath(args.save)), exist_ok=True)
        figs[0].savefig(args.save, dpi=120)
        print("saved (fallback) ->", args.save)

    if args.no_show:
        print("--no-show: figures computed, closing without displaying.")
        for f in figs:
            plt.close(f)
        return

    print("\nInteractive window(s) open — zoom/pan with the toolbar; "
          "close the window(s) to end.")
    plt.show(block=True)


if __name__ == "__main__":
    main()
