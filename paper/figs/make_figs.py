#!/usr/bin/env python3
"""Figures for the paper.

Reads only the shipped data in ../../certificate/ (no recomputation):
  implementation_A/exact_F_5e6.out            (1560 N, implementation A)
  implementation_B/exact_F_2e5_6e6.tsv        (409 N, implementation B)
  implementation_B/zeros_check_N25_58_74.out  (T_K vs K at N = 25, 58, 74)
Writes F1_fit.pdf, F2_residual.pdf, F3_zeros.pdf (vector) next to this script.
Usage:  python3 -I make_figs.py      (from any directory)
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CERT = HERE.parent.parent / "certificate"
C3 = CERT / "implementation_A/exact_F_5e6.out"
J3 = CERT / "implementation_B/exact_F_2e5_6e6.tsv"
J2 = CERT / "implementation_B/zeros_check_N25_58_74.out"

# Fixed categorical order (same series = same colour); references in grey.
BLUE, ORANGE, GREEN, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
GREY, INK, MUTED = "#52514e", "#000000", "#52514e"
CM = 1 / 2.54
W_FULL, W_085, H = 16.0 * CM, 13.6 * CM, 6.0 * CM
C0 = -4 * math.log(2 * math.pi)            # predicted constant -4 log(2 pi)

plt.rcParams.update({
    "font.family": "STIXGeneral", "mathtext.fontset": "stix",
    "font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GREY, "axes.labelcolor": INK,
    "xtick.color": GREY, "ytick.color": GREY,
    "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.5,
    "axes.axisbelow": True, "lines.linewidth": 1.75,
    "legend.frameon": False, "pdf.fonttype": 42,
})


def load_c3(path):
    """Columns: N F/N P^2/N (F-P^2)/N T/N rG(2N)/(2N) special."""
    a = np.loadtxt(path, comments="#")
    return a[:, 0], a[:, 1], a[:, 2]


def load_j3(path):
    """Columns: N F_over_N P2_over_N T_over_N rG2N_over_2N P (header line)."""
    a = np.loadtxt(path, skiprows=1)
    return a[:, 0], a[:, 1], a[:, 2]


def load_j2(path):
    """Parse 'N=.. T_stieltjes=..' blocks and their 'K=.. diff_vs_stieltjes=..' lines."""
    out, cur = {}, None
    for line in path.read_text().splitlines():
        m = re.match(r"N=(\d+) T_stieltjes=(\S+)", line)
        if m:
            cur = int(m.group(1))
            out[cur] = {"T": float(m.group(2)), "K": [], "TK": []}
            continue
        m = re.search(r"K=(\d+) .*T_K=(\S+)", line)
        if m and cur is not None:
            out[cur]["K"].append(int(m.group(1)))
            out[cur]["TK"].append(float(m.group(2)))
    return out


def fit(x, y):
    beta, c = np.polyfit(x, y, 1)
    return c, beta, float(np.std(y - (c + beta * x)))


def fig1(c3, j3):
    fig, ax = plt.subplots(figsize=(W_FULL, H))
    ax.scatter(c3[2], c3[1], s=6, color=BLUE, alpha=0.55, lw=0,
               label="implementation A (1560 $N$)")
    ax.scatter(j3[2], j3[1], s=10, color=ORANGE, alpha=0.8, lw=0, marker="D",
               label="implementation B (409 $N$)")
    xs = np.linspace(0, max(c3[2].max(), j3[2].max()) * 1.04, 50)
    ax.plot(xs, C0 + xs, color=GREY, lw=1.5, ls="--",
            label=r"$-4\log 2\pi + P(N)^2/N$")
    ax.set_xlabel(r"$P(N)^2/N$")
    ax.set_ylabel(r"$F(N)/N$")
    ax.set_xlim(0, xs[-1])
    ax.legend(loc="upper left", markerscale=2.0)
    ax.text(xs[-1] * 0.985, C0 + xs[-1] * 0.985 - 0.03, "slope 1", color=GREY,
            ha="right", va="top", fontsize=8)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "F1_fit.pdf")
    plt.close(fig)


def fig2(c3, j3):
    fig, ax = plt.subplots(figsize=(W_FULL, H))
    ax.scatter(c3[0], c3[1] - c3[2], s=6, color=BLUE, alpha=0.55, lw=0,
               label="implementation A")
    ax.scatter(j3[0], j3[1] - j3[2], s=10, color=ORANGE, alpha=0.8, lw=0,
               marker="D", label="implementation B")
    ax.axhline(C0, color=GREY, lw=1.5, ls="--", label=r"$-4\log 2\pi$")
    ax.set_xlabel(r"$N$ ($\times 10^{6}$)")
    ax.set_ylabel(r"$(F(N)-P(N)^2)/N$")
    ax.set_ylim(C0 - 6e-4, C0 + 6e-4)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:g}"))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:.4f}".replace("-", "\u2212")))
    ax.legend(loc="upper right", ncol=3, markerscale=2.0)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "F2_residual.pdf")
    plt.close(fig)


def fig3(j2):
    fig, ax = plt.subplots(figsize=(W_085, H))
    ax.axhline(0, color=GREY, lw=1.5, ls="--")
    lam25 = -math.log(5) ** 2 / 4
    ax.axhline(lam25, color=GREY, lw=1.0, ls=":")
    ax.text(7.6e2, lam25 - 0.03, r"$-\Lambda(25)^2/4$", color=GREY,
            ha="left", va="top", fontsize=8)
    dy = {58: 7, 74: -7, 25: -7}
    for n, col, mk in ((58, BLUE, "o"), (74, ORANGE, "D"), (25, GREEN, "s")):
        d = j2[n]
        diff = np.array(d["TK"]) - d["T"]
        ax.plot(d["K"], diff, color=col, marker=mk, ms=5, lw=1.75,
                label=f"$N={n}$")
        ax.annotate(f"$N={n}$", (d["K"][-1], diff[-1]), xytext=(7, dy[n]),
                    textcoords="offset points", va="center", fontsize=8,
                    color=INK)
    ax.set_xscale("log")
    ax.set_xlim(7e2, 1.0e5)
    ax.set_xlabel(r"number of zeros $K$")
    ax.set_ylabel(r"$T_K(N)-T(N)$")
    ax.legend(loc="lower right", ncol=3)
    fig.tight_layout(pad=0.3)
    fig.savefig(HERE / "F3_zeros.pdf")
    plt.close(fig)


def main():
    for p in (C3, J3, J2):
        if not p.is_file():
            sys.exit(f"missing input: {p}")
    c3, j3, j2 = load_c3(C3), load_j3(J3), load_j2(J2)
    for name, d in (("A (c3)", c3), ("B (j3)", j3)):
        c, beta, sd = fit(d[2], d[1])
        r = d[1] - d[2]
        print(f"{name}: n={len(d[0])} c={c:.6f} beta={beta:.6f} sd={sd:.2e} "
              f"mean(F-P^2)/N={r.mean():.6f} sd={r.std():.2e}")
    print(f"-4 log 2pi = {C0:.6f}")
    for n in (25, 58, 74):
        print(f"N={n} T={j2[n]['T']} K={j2[n]['K']} T_K={j2[n]['TK']}")
    fig1(c3, j3)
    fig2(c3, j3)
    fig3(j2)
    print("wrote F1_fit.pdf F2_residual.pdf F3_zeros.pdf")


if __name__ == "__main__":
    main()
