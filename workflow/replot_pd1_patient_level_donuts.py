#!/usr/bin/env python3
"""Create patient-level PD1 composition donuts.

Patient IDs follow the restored ACED convention: the first two numeric
filename tokens (for example, 016_1) identify a patient. Regions with the
same patient ID are pooled within patient before calculating one fraction.
The displayed percentage is the equal-weight mean of those patient fractions.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

ROOT = Path("significance_axons_immune_cells/paper_ready_tcell_marker_composition")
REGION_TABLE = ROOT / "marker_positive_fraction_region_level.csv"
OUTPUT_DIR = Path("manuscript final figures/figures")
OUTPUT_BASE = OUTPUT_DIR / "PD1_patient_level_composition_donuts"
SUMMARY_PATH = Path("manuscript final figures/tables/PD1_patient_level_composition_summary.csv")
PATIENT_PATH = Path("manuscript final figures/tables/PD1_patient_level_fractions.csv")

MPLCONFIG_DIR = ROOT / "mplconfig"
MPLCONFIG_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


MARKER = "PD1"
RADIUS_ORDER = [20, 50, 100]
ROW_ORDER = [
    ("CXCR4+ / CD44+ / FYN+ Axon", "NonProgressor"),
    ("CXCR4+ / CD44+ / FYN+ Axon", "Progressor"),
    ("NRXN1+ / MAOA+ Axon", "NonProgressor"),
    ("NRXN1+ / MAOA+ Axon", "Progressor"),
]
ROW_LABELS = {
    ("CXCR4+ / CD44+ / FYN+ Axon", "NonProgressor"): "CXCR4+ CD44+ FYN+ axon\nNon-progressor",
    ("CXCR4+ / CD44+ / FYN+ Axon", "Progressor"): "CXCR4+ CD44+ FYN+ axon\nProgressor",
    ("NRXN1+ / MAOA+ Axon", "NonProgressor"): "NRXN1+ MAOA+ axon\nNon-progressor",
    ("NRXN1+ / MAOA+ Axon", "Progressor"): "NRXN1+ MAOA+ axon\nProgressor",
}
POS_COLOR = "#C04A69"
NEG_COLOR = "#D8DCE1"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "DejaVu Sans",
        }
    )


def save_all(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


from patient_ids import patient_id_from_region


def main() -> None:
    setup_style()
    region = pd.read_csv(REGION_TABLE)
    data = region.loc[
        region["marker"].eq(MARKER)
        & region["progression"].isin(["NonProgressor", "Progressor"])
        & region["radius_um"].isin(RADIUS_ORDER)
    ].copy()
    data["patient_id"] = data["unique_region"].map(patient_id_from_region)

    # One observation per patient, axon population and radius. Regions from a
    # patient are pooled before calculating that patient's composition.
    patient = (
        data.groupby(["patient_id", "progression", "axon", "radius_um"], observed=True)
        .agg(
            n_regions=("unique_region", "nunique"),
            positive=("raw_positive_near", "sum"),
            negative=("raw_negative_near", "sum"),
            total=("raw_marker_split_near_total", "sum"),
        )
        .reset_index()
    )
    patient["pct_marker_positive_patient"] = 100.0 * patient["positive"] / patient["total"]
    if patient["pct_marker_positive_patient"].isna().any():
        raise ValueError("A patient-level composition has no classified PD1 observations.")

    summary = (
        patient.groupby(["axon", "progression", "radius_um"], observed=True)
        .agg(
            n_patients=("patient_id", "nunique"),
            n_regions=("n_regions", "sum"),
            mean_pct_patient=("pct_marker_positive_patient", "mean"),
            sd_pct_patient=("pct_marker_positive_patient", "std"),
            min_pct_patient=("pct_marker_positive_patient", "min"),
            max_pct_patient=("pct_marker_positive_patient", "max"),
        )
        .reset_index()
    )
    expected = {
        ("CXCR4+ / CD44+ / FYN+ Axon", "NonProgressor"): 7,
        ("CXCR4+ / CD44+ / FYN+ Axon", "Progressor"): 5,
        ("NRXN1+ / MAOA+ Axon", "NonProgressor"): 7,
    }
    for key, n_expected in expected.items():
        got = summary.loc[
            summary["axon"].eq(key[0]) & summary["progression"].eq(key[1]), "n_patients"
        ]
        if not got.eq(n_expected).all():
            raise ValueError(f"Unexpected patient count for {key}: {got.tolist()}")
    nrxn_progressor = summary.loc[
        summary["axon"].eq("NRXN1+ / MAOA+ Axon") & summary["progression"].eq("Progressor")
    ]
    if set(nrxn_progressor["n_patients"]) != {4, 5}:
        raise ValueError(f"Unexpected NRXN1+/MAOA+ progressor patient counts: {nrxn_progressor['n_patients'].tolist()}")
    # NRXN1+/MAOA+ at 20 um has no classified nearby T cells for patient 019_1,
    # so its progressor count is four rather than five.
    nrxn20 = summary.loc[
        summary["axon"].eq("NRXN1+ / MAOA+ Axon") & summary["radius_um"].eq(20)
    ]
    if set(nrxn20["n_patients"]) != {7, 4}:
        raise ValueError("Unexpected NRXN1+/MAOA+ 20-um patient counts")

    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SUMMARY_PATH, index=False)
    patient.to_csv(PATIENT_PATH, index=False)

    fig = plt.figure(figsize=(16.4, 18.2))
    gs = fig.add_gridspec(
        len(ROW_ORDER),
        len(RADIUS_ORDER) + 1,
        width_ratios=[1.50, 1.0, 1.0, 1.0],
        hspace=0.24,
        wspace=0.36,
    )
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=17,
               markerfacecolor=POS_COLOR, markeredgecolor="black", markeredgewidth=1.8,
               label="PD1+ T cells"),
        Line2D([0], [0], marker="o", linestyle="", markersize=17,
               markerfacecolor=NEG_COLOR, markeredgecolor="black", markeredgewidth=1.8,
               label="PD1- T cells"),
    ]
    fig.legend(handles=legend_handles, frameon=False, loc="upper center",
               bbox_to_anchor=(0.59, 0.975), ncol=2, prop={"size": 28},
               columnspacing=2.6, handletextpad=0.9)
    for x_position, radius in zip([0.455, 0.682, 0.91], RADIUS_ORDER):
        fig.text(x_position, 0.905, f"{radius} um", ha="center", va="bottom", fontsize=32)

    for row_index, (axon, progression) in enumerate(ROW_ORDER):
        label_ax = fig.add_subplot(gs[row_index, 0])
        label_ax.axis("off")
        label_ax.text(0.98, 0.5, ROW_LABELS[(axon, progression)], ha="right", va="center",
                      fontsize=26, linespacing=1.1)
        for col_index, radius in enumerate(RADIUS_ORDER, start=1):
            ax = fig.add_subplot(gs[row_index, col_index])
            ax.set_aspect("equal")
            ax.axis("off")
            row = summary.loc[
                summary["axon"].eq(axon)
                & summary["progression"].eq(progression)
                & summary["radius_um"].eq(radius)
            ]
            if row.empty:
                continue
            values = row.iloc[0]
            pct = float(values["mean_pct_patient"])
            ax.pie([pct, 100.0 - pct], colors=[POS_COLOR, NEG_COLOR], startangle=90,
                   counterclock=False, radius=1.0,
                   wedgeprops={"width": 0.38, "edgecolor": "white", "linewidth": 2.8})
            ax.text(0, 0, f"{int(round(pct))}%\nPD1+", ha="center", va="center",
                    fontsize=25, linespacing=1.0)
            ax.text(0.5, -0.15,
                    f"n={int(values['n_patients'])} patients\n"
                    f"n={int(values['n_regions'])} lesions",
                    transform=ax.transAxes, ha="center", va="top", fontsize=16, color="#555555")

    fig.text(
        0.5,
        0.018,
        "Each patient contributes one fraction; multiple lesions are pooled within patient.",
        ha="center",
        va="bottom",
        fontsize=14,
        color="#555555",
    )
    fig.subplots_adjust(left=0.035, right=0.99, bottom=0.07, top=0.85)
    save_all(fig, OUTPUT_BASE)
    print(OUTPUT_BASE.with_suffix(".pdf"))
    print(OUTPUT_BASE.with_suffix(".png"))
    print(OUTPUT_BASE.with_suffix(".svg"))
    print(SUMMARY_PATH)
    print(PATIENT_PATH)


if __name__ == "__main__":
    main()
