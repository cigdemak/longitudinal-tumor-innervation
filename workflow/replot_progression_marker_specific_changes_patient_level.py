#!/usr/bin/env python3
"""Make patient-level versions of the progression marker-effect panels.

The original panels use one row per lesion/region.  This sensitivity replot
uses the ACED filename convention used for this reproduction: the first two
numeric filename tokens identify the patient (for example, ``016_1``).  All
lesions with the same patient ID are pooled before one value is calculated for
that patient.  The original region-level figures are not overwritten.

The abundance panel is a descriptive patient-level re-expression of the raw
nearby-count endpoint.  For each patient it pools marker-positive nearby-cell
counts and tissue totals, then reports nearby marker-positive cells per 1,000
tissue marker-positive cells.  The plotted effect is the difference between
the Progressor and NonProgressor means.  The PD1 panel reports the difference
in the equal-weight patient means of the pooled PD1-positive fractions, in
percentage points.  Mann-Whitney p values and FDR q values are recalculated at
the patient level and are saved with the source values.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path.cwd()
ABUNDANCE_SOURCE = ROOT / "significance_axons_immune_cells/paper_ready_tcell_marker_positivity/tcell_marker_region_level_metrics.csv"
COMPOSITION_SOURCE = ROOT / "significance_axons_immune_cells/paper_ready_tcell_marker_composition/marker_positive_fraction_region_level.csv"
FIGURE_DIR = ROOT / "manuscript final figures/figures"
TABLE_DIR = ROOT / "manuscript final figures/tables"
LOG_DIR = ROOT / "manuscript final figures/logs/mplconfig"
LOG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(LOG_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


RADIUS_ORDER = [20, 50, 100]
PROGRESSION_ORDER = ["NonProgressor", "Progressor"]
ABUNDANCE_MARKERS = ["AR", "Ki67", "GZMB"]
ABUNDANCE_LABELS = {"AR": "AR+ T cells", "Ki67": "Ki67+ T cells", "GZMB": "GZMB+ T cells"}
ABUNDANCE_COLORS = {"AR": "#287bb7", "Ki67": "#ff7f0e", "GZMB": "#2ca02c"}
AXON_ORDER = ["CXCR4+ / CD44+ / FYN+ Axon", "NRXN1+ / MAOA+ Axon"]
AXON_LABELS = {
    "CXCR4+ / CD44+ / FYN+ Axon": "CXCR4+ / CD44+ / FYN+",
    "NRXN1+ / MAOA+ Axon": "NRXN1+ / MAOA+",
}
AXON_COLORS = {
    "CXCR4+ / CD44+ / FYN+ Axon": "#7b3294",
    "NRXN1+ / MAOA+ Axon": "#c51b7d",
}

PATIENT_RE = re.compile(r"^/Batch\d+/([^_]+_[^_]+)_")


from patient_ids import patient_id_from_region


def setup_style() -> None:
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
            "axes.linewidth": 2.2,
            "xtick.major.width": 2.2,
            "ytick.major.width": 2.2,
            "xtick.major.size": 8,
            "ytick.major.size": 8,
            "lines.solid_capstyle": "round",
        }
    )


def save_figure(fig: plt.Figure, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def q_stars(q_value: object) -> str:
    q = pd.to_numeric(pd.Series([q_value]), errors="coerce").iloc[0]
    if not np.isfinite(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def patient_level_tests(
    patient: pd.DataFrame,
    value_col: str,
    group_cols: list[str],
    effect_name: str,
    scale: str,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for keys, group in patient.groupby(group_cols, observed=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        groups = {
            status: pd.to_numeric(group.loc[group["progression"].eq(status), value_col], errors="coerce").dropna()
            for status in PROGRESSION_ORDER
        }
        p_values = groups["Progressor"]
        n_values = groups["NonProgressor"]
        if len(p_values) and len(n_values):
            _, p_value = mannwhitneyu(p_values, n_values, alternative="two-sided")
        else:
            p_value = np.nan
        row = dict(zip(group_cols, keys))
        row.update(
            {
                "n_progressor_patients": int(len(p_values)),
                "n_nonprogressor_patients": int(len(n_values)),
                "n_progressor_lesions": int(
                    group.loc[group["progression"].eq("Progressor"), "n_lesions"].sum()
                ),
                "n_nonprogressor_lesions": int(
                    group.loc[group["progression"].eq("NonProgressor"), "n_lesions"].sum()
                ),
                "mean_progressor": float(p_values.mean()) if len(p_values) else np.nan,
                "mean_nonprogressor": float(n_values.mean()) if len(n_values) else np.nan,
                "median_progressor": float(p_values.median()) if len(p_values) else np.nan,
                "median_nonprogressor": float(n_values.median()) if len(n_values) else np.nan,
                "effect_progressor_minus_nonprogressor": (
                    float(p_values.mean() - n_values.mean()) if len(p_values) and len(n_values) else np.nan
                ),
                "p_value_two_sided": float(p_value) if np.isfinite(p_value) else np.nan,
                "effect_name": effect_name,
                "effect_scale": scale,
            }
        )
        rows.append(row)
    result = pd.DataFrame(rows)
    if not result.empty:
        valid = result["p_value_two_sided"].notna()
        result["q_value_fdr"] = np.nan
        if valid.any():
            result.loc[valid, "q_value_fdr"] = multipletests(
                result.loc[valid, "p_value_two_sided"], method="fdr_bh"
            )[1]
    return result


def aggregate_abundance() -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(ABUNDANCE_SOURCE)
    data = data.loc[
        data["df"].eq("aced")
        & data["marker_state"].eq("positive")
        & data["axon"].eq("NRXN1+ / MAOA+ Axon")
        & data["marker"].isin(ABUNDANCE_MARKERS)
        & data["progression"].isin(PROGRESSION_ORDER)
        & data["radius_um"].isin(RADIUS_ORDER)
    ].copy()
    data["patient_id"] = data["unique_region"].map(patient_id_from_region)
    patient = (
        data.groupby(["patient_id", "progression", "marker", "radius_um"], observed=True)
        .agg(
            n_lesions=("unique_region", "nunique"),
            nearby_marker_positive=("raw_n_near", "sum"),
            tissue_total=("tissue_total", "sum"),
        )
        .reset_index()
    )
    patient["nearby_per_1000_tissue_marker_positive"] = np.where(
        patient["tissue_total"] > 0,
        1000.0 * patient["nearby_marker_positive"] / patient["tissue_total"],
        np.nan,
    )
    tests = patient_level_tests(
        patient,
        "nearby_per_1000_tissue_marker_positive",
        ["marker", "radius_um"],
        "mean Progressor minus NonProgressor",
        "nearby marker-positive T cells per 1,000 tissue marker-positive T cells",
    )
    tests["axon"] = "NRXN1+ / MAOA+ Axon"
    return patient, tests


def aggregate_pd1_fraction() -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(COMPOSITION_SOURCE)
    data = data.loc[
        data["df"].eq("aced")
        & data["marker"].eq("PD1")
        & data["axon"].isin(AXON_ORDER)
        & data["progression"].isin(PROGRESSION_ORDER)
        & data["radius_um"].isin(RADIUS_ORDER)
    ].copy()
    data["patient_id"] = data["unique_region"].map(patient_id_from_region)
    patient = (
        data.groupby(["patient_id", "progression", "axon", "radius_um"], observed=True)
        .agg(
            n_lesions=("unique_region", "nunique"),
            pd1_positive=("raw_positive_near", "sum"),
            pd1_negative=("raw_negative_near", "sum"),
            nearby_t_cells=("raw_marker_split_near_total", "sum"),
        )
        .reset_index()
    )
    patient["pd1_positive_fraction_pct"] = np.where(
        patient["nearby_t_cells"] > 0,
        100.0 * patient["pd1_positive"] / patient["nearby_t_cells"],
        np.nan,
    )
    tests = patient_level_tests(
        patient,
        "pd1_positive_fraction_pct",
        ["axon", "radius_um"],
        "mean Progressor minus NonProgressor",
        "percentage-point difference in equal-weight patient PD1-positive fractions",
    )
    return patient, tests


def finish_axis(ax: plt.Axes, y_values: pd.Series, ylabel: str, star_pad: float) -> None:
    ax.axhline(0, color="#666666", linewidth=2.0, zorder=1)
    ax.set_xlim(14, 106)
    ax.set_xticks(RADIUS_ORDER)
    ax.set_xlabel("Radius (um)", fontsize=38, labelpad=12)
    ax.set_ylabel(ylabel, fontsize=34, labelpad=12)
    ax.tick_params(axis="both", labelsize=34)
    ax.grid(axis="y", color="#d8d8d8", linewidth=1.1)
    finite = pd.to_numeric(y_values, errors="coerce").dropna().to_numpy(dtype=float)
    if finite.size:
        lower = min(float(finite.min()), 0.0)
        upper = max(float(finite.max()), 0.0)
        span = max(upper - lower, abs(star_pad) * 4.0, 1.0)
        ax.set_ylim(lower - 0.12 * span, upper + 0.14 * span)
    else:
        ax.set_ylim(-1, 1)
    sns_despine(ax)


def sns_despine(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def line_legend(fig: plt.Figure, labels: list[str], colors: dict[str, str], ncol: int = 1) -> None:
    handles = [
        Line2D([0], [0], color=colors[label], marker="o", markersize=12, linewidth=4.5, label=label)
        for label in labels
    ]
    fig.legend(
        handles=handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.12),
        ncol=ncol,
        prop={"size": 26},
        handlelength=2.0,
        handletextpad=0.7,
        columnspacing=1.2,
    )


def plot_abundance(patient: pd.DataFrame, tests: pd.DataFrame, base: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 9.2))
    for marker in ABUNDANCE_MARKERS:
        grp = tests.loc[tests["marker"].eq(marker)].sort_values("radius_um")
        if grp.empty:
            continue
        ax.plot(
            grp["radius_um"],
            grp["effect_progressor_minus_nonprogressor"],
            marker="o",
            markersize=13,
            linewidth=4.5,
            color=ABUNDANCE_COLORS[marker],
            label=ABUNDANCE_LABELS[marker],
            zorder=3,
        )
        for _, row in grp.iterrows():
            star = q_stars(row["q_value_fdr"])
            if star:
                ax.text(
                    row["radius_um"],
                    row["effect_progressor_minus_nonprogressor"] + 0.05 * max(abs(tests["effect_progressor_minus_nonprogressor"].min()), 1),
                    star,
                    ha="center",
                    va="bottom",
                    fontsize=34,
                    color=ABUNDANCE_COLORS[marker],
                )
    finish_axis(
        ax,
        tests["effect_progressor_minus_nonprogressor"],
        "Progressor effect\n(per 1,000 cells)",
        2.0,
    )
    line_legend(fig, [ABUNDANCE_LABELS[m] for m in ABUNDANCE_MARKERS], {ABUNDANCE_LABELS[m]: ABUNDANCE_COLORS[m] for m in ABUNDANCE_MARKERS})
    fig.text(
        0.5,
        0.205,
        "n=7 NonProgressor / 5 Progressor patients; 15 / 8 lesions at each radius",
        ha="center",
        va="bottom",
        fontsize=12,
        color="#555555",
    )
    fig.tight_layout(rect=[0.13, 0.22, 1, 1])
    save_figure(fig, base)


def plot_pd1(patient: pd.DataFrame, tests: pd.DataFrame, base: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 9.2))
    for axon in AXON_ORDER:
        grp = tests.loc[tests["axon"].eq(axon)].sort_values("radius_um")
        if grp.empty:
            continue
        ax.plot(
            grp["radius_um"],
            grp["effect_progressor_minus_nonprogressor"],
            marker="o",
            markersize=13,
            linewidth=4.5,
            color=AXON_COLORS[axon],
            label=AXON_LABELS[axon],
            zorder=3,
        )
        for _, row in grp.iterrows():
            star = q_stars(row["q_value_fdr"])
            if star:
                ax.text(
                    row["radius_um"],
                    row["effect_progressor_minus_nonprogressor"] + 1.6,
                    star,
                    ha="center",
                    va="bottom",
                    fontsize=34,
                    color=AXON_COLORS[axon],
                )
    finish_axis(
        ax,
        tests["effect_progressor_minus_nonprogressor"],
        "Progressor effect\n(percentage points)",
        1.6,
    )
    line_legend(fig, [AXON_LABELS[a] for a in AXON_ORDER], {AXON_LABELS[a]: AXON_COLORS[a] for a in AXON_ORDER})
    fig.text(
        0.5,
        0.205,
        "CXCR4: 7/5 patients, 15/8 lesions; NRXN1: 20 um 7/4, 11/7; 50 um 7/5, 12/8; 100 um 7/5, 13/8",
        ha="center",
        va="bottom",
        fontsize=11.5,
        color="#555555",
    )
    fig.tight_layout(rect=[0.13, 0.22, 1, 1])
    save_figure(fig, base)


def main() -> None:
    setup_style()
    abundance_patient, abundance_tests = aggregate_abundance()
    pd1_patient, pd1_tests = aggregate_pd1_fraction()

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    abundance_patient.to_csv(TABLE_DIR / "progression_marker_patient_level_abundance.csv", index=False)
    abundance_tests.to_csv(TABLE_DIR / "progression_marker_patient_level_abundance_tests.csv", index=False)
    pd1_patient.to_csv(TABLE_DIR / "progression_marker_patient_level_pd1_fraction.csv", index=False)
    pd1_tests.to_csv(TABLE_DIR / "progression_marker_patient_level_pd1_fraction_tests.csv", index=False)

    abundance_base = FIGURE_DIR / "fig_progression_marker_specific_changes_abundance_patient_level"
    pd1_base = FIGURE_DIR / "fig_progression_marker_specific_changes_pd1_fraction_patient_level"
    plot_abundance(abundance_patient, abundance_tests, abundance_base)
    plot_pd1(pd1_patient, pd1_tests, pd1_base)

    for path in [
        abundance_base.with_suffix(".pdf"),
        abundance_base.with_suffix(".png"),
        abundance_base.with_suffix(".svg"),
        pd1_base.with_suffix(".pdf"),
        pd1_base.with_suffix(".png"),
        pd1_base.with_suffix(".svg"),
    ]:
        print(path)


if __name__ == "__main__":
    main()
