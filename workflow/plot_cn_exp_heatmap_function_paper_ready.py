#!/usr/bin/env python3
"""Plot CN cell-type heatmaps with the provided cn_exp_heatmap logic."""

from __future__ import annotations

import argparse
import os
import pathlib
from pathlib import Path


DEFAULT_H5ADS = [
    Path("data/adata_aced_immune_all_CN_k50_n10_with_epithelia_filtered.h5ad"),
    Path("data/adata_BLimmune_all_CN_k50_n10_with_epithelia.h5ad"),
]
DEFAULT_OUTDIR = Path("manuscript final figures/figures")
DEFAULT_TABLE_DIR = Path("manuscript final figures/tables")
DEFAULT_PREFIX = "fig_cn_exp_heatmap_function"

MPLCONFIG_DIR = Path("manuscript final figures/mplconfig")
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_hex


CN_ANNOTATION_COLORS = {
    0: "#7B8C63",
    1: "#275BC9",
    2: "#5EC165",
    3: "#F2A1A2",
    4: "#DA34A2",
    5: "#1C4879",
    6: "#F3AB3C",
    7: "#8E5790",
    8: "#EE662C",
    9: "#357F7F",
}

CELL_LABELS = {
    "CXCR4+ / CD44+ / FYN+ Axon": "CXCR4+ / CD44+ / FYN+ Axon",
    "NRXN1+ / MAOA+ Axon": "NRXN1+ / MAOA+ Axon",
    "TH+ Axon": "TH+ Axon",
    "CD90+ Axon": "CD90+ Axon",
    "CD90+ / INA+ Axon": "CD90+ / INA+ Axon",
    "T Cells": "T Cells",
    "CD4+ T Cells": "CD4+ T Cells",
    "FOXP3+ T Cells": "FOXP3+ T Cells",
    "Macrophages": "Macrophages",
    "CD20+": "CD20+",
    "CD66b+": "CD66b+",
    "Tryp+": "Tryp+",
    "AR+ Epithelia": "AR+ Epithelia",
    "Basal Epithelia": "Basal Epithelia",
    "AMACR-high Epithelia": "AMACR-high Epithelia",
    "aSMA+": "aSMA+",
    "AR+aSMA+": "AR+ aSMA+",
    "Other Stroma": "Other Stroma",
    "Mesenchymal/Endothelial": "Mesenchymal/Endothelial",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create paper-ready CN heatmaps with the provided cn_exp_heatmap function."
    )
    parser.add_argument("--h5ad", type=Path, nargs="+", default=DEFAULT_H5ADS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--cn-column", default="CN_k50_n10")
    parser.add_argument("--cell-type-column", default="cell_type_sync")
    parser.add_argument(
        "--cohort",
        choices=["all", "aced", "bl", "combined"],
        default="all",
        help="Which cohort to plot. 'all' writes ACED, BL, and combined.",
    )
    parser.add_argument("--no-row-cluster", action="store_true")
    parser.add_argument("--no-col-cluster", action="store_true")
    parser.add_argument(
        "--show-dendrograms",
        action="store_true",
        help="Show row and column dendrograms. Hidden by default for manuscript output.",
    )
    parser.add_argument("--figsize", type=float, nargs=2, default=(17.57, 12.5))
    parser.add_argument("--font-scale", type=float, default=1.5)
    return parser.parse_args()


def setup_style() -> None:
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.linewidth": 2.2,
            "xtick.major.width": 2.2,
            "ytick.major.width": 2.2,
        }
    )


def hf_generate_random_colors(n_colors: int, rand_seed: int = 1) -> list[str]:
    rng = np.random.default_rng(rand_seed)
    hues = np.linspace(0, 1, n_colors, endpoint=False)
    rng.shuffle(hues)
    colors = sns.husl_palette(n_colors, h=0.01, s=0.78, l=0.58)
    return [to_hex(color) for _, color in sorted(zip(hues, colors), key=lambda item: item[0])]


def infer_cohort(path: Path) -> str:
    name = path.name.lower()
    if "aced" in name:
        return "aced"
    if "bl" in name:
        return "bl"
    return path.stem


def clean_category_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for column in out.select_dtypes(["category"]).columns:
        out[column] = out[column].astype(str)
    return out


def cn_sort_key(label: object) -> tuple[int, object]:
    text = str(label).replace("CN", "").strip()
    try:
        return (0, int(float(text)))
    except ValueError:
        return (1, str(label))


def cn_number(label: object) -> int | None:
    text = str(label).replace("CN", "").strip()
    try:
        return int(float(text))
    except ValueError:
        return None


def fixed_cn_palette(cn_values: list[object], rand_seed: int = 1) -> dict[object, str]:
    ordered = sorted(cn_values, key=cn_sort_key)
    fallback = hf_generate_random_colors(len(ordered), rand_seed=rand_seed)
    palette = {}
    for i, value in enumerate(ordered):
        number = cn_number(value)
        palette[value] = CN_ANNOTATION_COLORS.get(number, fallback[i])
    return palette


def display_cell_type(label: object) -> str:
    return CELL_LABELS.get(str(label), str(label).replace("_", " "))


def read_obs(paths: list[Path], cn_column: str, cell_type_column: str) -> pd.DataFrame:
    frames = []
    for path in paths:
        cohort = infer_cohort(path)
        adata = ad.read_h5ad(path, backed="r")
        try:
            missing = [column for column in [cn_column, cell_type_column] if column not in adata.obs.columns]
            if missing:
                raise KeyError(f"{path} is missing required obs columns: {missing}")
            obs = clean_category_frame(adata.obs[[cn_column, cell_type_column]].copy())
        finally:
            adata.file.close()

        obs = obs.rename(columns={cn_column: "cellular_neighborhood", cell_type_column: "cell_type"})
        obs["cohort"] = cohort
        obs["cellular_neighborhood"] = obs["cellular_neighborhood"].astype(str).map(lambda value: f"CN {value}")
        obs["cell_type"] = obs["cell_type"].map(display_cell_type)
        obs = obs.loc[
            obs["cellular_neighborhood"].notna()
            & obs["cell_type"].notna()
            & ~obs["cellular_neighborhood"].astype(str).str.lower().isin(["", "nan", "none", "na"])
            & ~obs["cell_type"].astype(str).str.lower().isin(["", "nan", "none", "na"])
        ].copy()
        frames.append(obs[["cohort", "cellular_neighborhood", "cell_type"]])
    return pd.concat(frames, ignore_index=True)


def make_adata(obs: pd.DataFrame) -> ad.AnnData:
    obs = obs[["cellular_neighborhood", "cell_type"]].reset_index(drop=True).copy()
    obs.index = obs.index.astype(str)
    return ad.AnnData(obs=obs)


def cn_exp_heatmap(
    adata,
    cluster_col,
    cn_col,
    palette=None,
    savefig=False,
    output_fname="",
    output_dir="./",
    row_clus=True,
    col_clus=True,
    rand_seed=1,
    figsize=(10, 5),
    show_dendrograms=False,
    font_scale=1.0,
):
    """
    Create a heatmap of expression data, clustered by rows and columns.

    This is the provided CN heatmap calculation with manuscript styling added
    around the plot export. It returns the Seaborn ClusterGrid and the plotted
    log2 fold-change matrix for reproducibility.
    """
    data = adata.obs
    output_dir = pathlib.Path(output_dir)

    if palette is None:
        if cn_col + "_colors" not in adata.uns.keys():
            cn_values = sorted(adata.obs[cn_col].unique(), key=cn_sort_key)
            palette = fixed_cn_palette(cn_values, rand_seed=rand_seed)
            adata.uns[cn_col + "_colors"] = [palette[value] for value in cn_values]
        else:
            cn_values = sorted(adata.obs[cn_col].unique(), key=cn_sort_key)
            palette = dict(zip(cn_values, adata.uns[cn_col + "_colors"]))

    neigh_data = pd.DataFrame({cn_col: list(palette.keys()), "color": list(palette.values())})
    neigh_data.set_index(keys=cn_col, inplace=True)

    df3 = pd.concat([data, pd.get_dummies(data[cluster_col])], axis=1)
    sum_cols2 = pd.Index(df3[cluster_col].unique()).astype(str)
    values2 = df3[sum_cols2].values
    cell_list = sum_cols2.tolist()
    cell_list.append(cn_col)

    subset = df3[cell_list]
    niche_sub = subset.groupby(cn_col, observed=True).sum()
    niche_sub = niche_sub.reindex(sorted(niche_sub.index, key=cn_sort_key))
    niche_df = niche_sub.apply(lambda x: x / x.sum() * 10, axis=1)
    neigh_clusters = niche_df.to_numpy()

    tissue_avgs = values2.mean(axis=0)

    fc_2 = np.log2(
        ((neigh_clusters + tissue_avgs) / (neigh_clusters + tissue_avgs).sum(axis=1, keepdims=True))
        / tissue_avgs
    )
    fc_2 = pd.DataFrame(fc_2, columns=sum_cols2)
    fc_2.set_index(niche_df.index, inplace=True)

    s = sns.clustermap(
        fc_2,
        vmin=-3,
        vmax=3,
        cmap="vlag",
        figsize=figsize,
        row_cluster=row_clus,
        col_cluster=col_clus,
        row_colors=[neigh_data.reindex(fc_2.index)["color"]],
        cbar_pos=None,
        dendrogram_ratio=(0.1, 0.13),
        colors_ratio=0.03,
        linewidths=0,
        linecolor="white",
    )

    s.ax_row_dendrogram.set_visible(show_dendrograms and row_clus)
    s.ax_col_dendrogram.set_visible(show_dendrograms and col_clus)
    font_scale = float(font_scale)
    if not np.isfinite(font_scale) or font_scale <= 0:
        font_scale = 1.0
    x_tick_font = 20 * font_scale
    y_tick_font = 20 * font_scale

    s.ax_heatmap.set_xlabel("")
    s.ax_heatmap.set_ylabel("")
    s.ax_heatmap.tick_params(axis="x", labelsize=x_tick_font, length=0, pad=3 * font_scale)
    s.ax_heatmap.tick_params(axis="y", labelsize=y_tick_font, length=4, pad=20 * font_scale)
    s.ax_heatmap.set_xticklabels(
        s.ax_heatmap.get_xticklabels(),
        rotation=90,
        ha="center",
        fontweight="normal",
    )
    y_labels = []
    for tick in s.ax_heatmap.get_yticklabels():
        number = cn_number(tick.get_text())
        y_labels.append(str(number) if number is not None else tick.get_text())
    s.ax_heatmap.set_yticklabels(y_labels, rotation=0, fontweight="normal")
    for tick in s.ax_heatmap.get_xticklabels() + s.ax_heatmap.get_yticklabels():
        tick.set_fontweight("normal")
    s.ax_heatmap.yaxis.set_ticks_position("right")
    s.ax_heatmap.yaxis.set_label_position("left")

    if not show_dendrograms:
        s.ax_row_dendrogram.set_position((0.075, 0.610, 0.095, 0.330))
        s.ax_col_dendrogram.set_position((0.227, 0.942, 0.708, 0.001))
        s.ax_row_colors.set_position((0.200, 0.610, 0.027, 0.330))
        s.ax_heatmap.set_position((0.227, 0.610, 0.708, 0.330))

    if savefig:
        output_dir.mkdir(parents=True, exist_ok=True)
        for suffix in [".pdf", ".png", ".svg"]:
            kwargs = {}
            if suffix == ".png":
                kwargs["dpi"] = 300
            s.figure.savefig(output_dir / (output_fname + suffix), **kwargs)

    return s, fc_2


def cohorts_to_plot(obs: pd.DataFrame, choice: str) -> dict[str, pd.DataFrame]:
    if choice == "combined":
        return {"combined": obs.copy()}
    if choice == "all":
        out = {cohort: group.copy() for cohort, group in obs.groupby("cohort", observed=True)}
        out["combined"] = obs.copy()
        return out
    return {choice: obs.loc[obs["cohort"].eq(choice)].copy()}


def main() -> None:
    args = parse_args()
    setup_style()
    obs = read_obs(args.h5ad, args.cn_column, args.cell_type_column)
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    for cohort, cohort_obs in cohorts_to_plot(obs, args.cohort).items():
        if cohort_obs.empty:
            continue
        adata = make_adata(cohort_obs)
        output_fname = f"{args.prefix}_{cohort}"
        grid, matrix = cn_exp_heatmap(
            adata,
            cluster_col="cell_type",
            cn_col="cellular_neighborhood",
            savefig=True,
            output_fname=output_fname,
            output_dir=args.outdir,
            row_clus=not args.no_row_cluster,
            col_clus=not args.no_col_cluster,
            figsize=tuple(args.figsize),
            show_dendrograms=args.show_dendrograms,
            font_scale=args.font_scale,
        )
        matrix.to_csv(args.table_dir / f"{output_fname}_log2_fold_change_matrix.csv")
        plt.close(grid.figure)
        print(f"Wrote {args.outdir / (output_fname + '.pdf')}")


if __name__ == "__main__":
    main()
