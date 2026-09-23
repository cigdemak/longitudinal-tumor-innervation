#!/usr/bin/env python3
"""Create a standalone color-scale legend for the CN heatmap."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


DEFAULT_OUTDIR = Path("manuscript final figures/figures")
DEFAULT_PREFIX = "fig_cn_exp_heatmap_color_scale_legend"
MPLCONFIG_DIR = Path("manuscript final figures/mplconfig")
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a standalone CN heatmap colorbar legend.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--figsize", type=float, nargs=2, default=(8.8, 2.7))
    parser.add_argument("--font-size", type=float, default=36)
    parser.add_argument("--label-font-size", type=float, default=36)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--vmin", type=float, default=-3)
    parser.add_argument("--vmax", type=float, default=3)
    parser.add_argument("--ticks", type=float, nargs="+", default=[-3, 0, 3])
    return parser.parse_args()


def setup_style() -> None:
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
            "font.weight": "normal",
            "axes.linewidth": 2.6,
            "xtick.major.width": 2.6,
            "ytick.major.width": 2.6,
        }
    )


def main() -> None:
    args = parse_args()
    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=tuple(args.figsize))
    ax = fig.add_axes([0.08, 0.58, 0.84, 0.28])
    cmap = sns.color_palette("vlag", as_cmap=True)
    norm = Normalize(vmin=args.vmin, vmax=args.vmax)
    mapper = ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array([])

    cbar = fig.colorbar(mapper, cax=ax, orientation="horizontal", ticks=args.ticks)
    cbar.set_label("")
    cbar.ax.tick_params(
        axis="x",
        labelsize=args.font_size,
        width=2.6,
        length=10,
        pad=8,
    )
    cbar.outline.set_linewidth(2.6)

    for tick in cbar.ax.get_xticklabels():
        tick.set_fontweight("normal")

    fig.text(
        0.5,
        0.17,
        "log2 fold change\nover tissue average",
        ha="center",
        va="center",
        fontsize=args.label_font_size,
        linespacing=0.95,
    )

    for suffix in [".pdf", ".png", ".svg"]:
        kwargs = {}
        if suffix == ".png":
            kwargs["dpi"] = args.dpi
        fig.savefig(args.outdir / f"{args.prefix}{suffix}", **kwargs)
    plt.close(fig)
    print(f"Wrote {args.outdir / (args.prefix + '.pdf')}")


if __name__ == "__main__":
    main()
