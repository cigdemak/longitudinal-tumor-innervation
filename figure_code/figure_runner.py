from __future__ import annotations

import os
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
WORKFLOW = ROOT / "workflow"
H5AD_FILES = {
    "cn_aced": "adata_aced_immune_all_CN_k50_n10_with_epithelia_filtered.h5ad",
    "cn_bl": "adata_BLimmune_all_CN_k50_n10_with_epithelia.h5ad",
    "marker_aced": "adata_aced_immune_updated_tcell_marker_celltypes.h5ad",
    "marker_bl": "adata_bl_immune_updated_tcell_marker_celltypes.h5ad",
}
AXON_TAGS = (
    "CXCR4pos___CD44pos___FYNpos_Axon",
    "NRXN1pos___MAOApos_Axon",
)
PPA_FILES = {
    f"{marker}_{axon}": f"aced_immune_{marker}_{axon}_ppa.h5ad"
    for marker in ("GZMB", "PD1")
    for axon in AXON_TAGS
}
PANEL_FILES = {
    "figure_4a_cn_heatmap": ("cn_aced", "cn_bl"),
    "figure_4bc_whole_tissue": ("cn_aced", "cn_bl"),
    "figure_4e_window_density": ("cn_aced", "cn_bl"),
    "figure_5b_cn_progression": ("cn_aced",),
    "figure_5d_local_redistribution": ("local_progression_statistics.csv",),
    "figure_6a_pd1_composition": tuple(f"PD1_{axon}" for axon in AXON_TAGS),
    "figure_6b_whole_tissue_markers": ("marker_aced", "marker_bl"),
    "figure_6cd_axon_marker_fractions": tuple(PPA_FILES),
    "figure_6ef_progression_effects": ("local_progression_statistics.csv",),
}
FIGURE_OUTPUTS = {
    "figure_4a_cn_heatmap": ("fig_cn_cell_type_heatmap_combined.pdf",),
    "figure_4bc_whole_tissue": (
        "fig_whole_tissue_grade_composition_separate_panel.pdf",
        "fig_whole_tissue_grade_region_immune_fraction_separate_panel.pdf",
        "fig_whole_tissue_grade_immune_density_separate_panel.pdf",
    ),
    "figure_4e_window_density": ("fig_axon_immune_target_window_density_background_vs_axon_aced_bl_combined.pdf",),
    "figure_5b_cn_progression": ("aced_progression_paper_ready.pdf",),
    "figure_5d_local_redistribution": ("fig_progression_local_redistribution.pdf",),
    "figure_6a_pd1_composition": ("PD1_raw_composition_donuts.pdf",),
    "figure_6b_whole_tissue_markers": (
        "fig_whole_tissue_tcell_marker_grade_whole_sample_fraction_separate_panel.pdf",
        "fig_whole_tissue_tcell_marker_grade_region_fraction_separate_panel.pdf",
        "fig_whole_tissue_tcell_marker_grade_positive_tcell_density_separate_panel.pdf",
        "fig_whole_tissue_tcell_marker_progression_whole_sample_fraction_separate_panel.pdf",
        "fig_whole_tissue_tcell_marker_progression_region_fraction_separate_panel.pdf",
        "fig_whole_tissue_tcell_marker_progression_positive_tcell_density_separate_panel.pdf",
    ),
    "figure_6cd_axon_marker_fractions": (
        "PD1_positive_fraction_barplot.pdf",
        "GZMB_positive_fraction_barplot.pdf",
    ),
    "figure_6ef_progression_effects": (
        "fig_progression_marker_specific_changes_pd1_fraction_separate_panel.pdf",
        "fig_progression_marker_specific_changes_abundance_separate_panel.pdf",
    ),
}


def locate_input(data_dir: Path, name: str) -> Path | None:
    if name == "local_progression_statistics.csv":
        candidates = (
            data_dir / "focused_fdr_master_analysis_table.csv",
            data_dir / "archived" / "focused_fdr_master_analysis_table.csv",
            data_dir / "manuscript_panels_results" / "tables" / "focused_fdr_master_analysis_table.csv",
            data_dir.parent / "manuscript_panels_results" / "tables" / "focused_fdr_master_analysis_table.csv",
        )
    else:
        filename = {**H5AD_FILES, **PPA_FILES}.get(name, name)
        candidates = (
            data_dir / filename,
            data_dir / "data" / filename,
            data_dir / "PPA_axon_subtypes_marker_split_T_cells" / filename,
            data_dir.parent / "PPA_axon_subtypes_marker_split_T_cells" / filename,
        )
    return next((path for path in candidates if path.is_file()), None)


def panel_commands(panel: str, inputs: dict[str, Path], archive_target: Path) -> list[list[str]]:
    final = "manuscript final figures/figures"
    tables = "manuscript final figures/tables"
    grade_tables = "manuscript_panels_results/whole_tissue_grade/combined"
    marker_grade_tables = "manuscript_panels_results/whole_tissue_tcell_marker_grade/combined"
    marker_progression_tables = "manuscript_panels_results/whole_tissue_tcell_marker_progression"
    cn_tables = "results/cn_neighborhood_enrichment_progression_aced"
    cn_aced = str(inputs.get(H5AD_FILES["cn_aced"], ""))
    cn_bl = str(inputs.get(H5AD_FILES["cn_bl"], ""))
    marker_aced = str(inputs.get(H5AD_FILES["marker_aced"], ""))
    marker_bl = str(inputs.get(H5AD_FILES["marker_bl"], ""))
    ppa_files = [inputs[name] for name in PPA_FILES.values() if name in inputs]
    ppa_dir = str(ppa_files[0].parent) if ppa_files else ""
    if panel == "figure_4a_cn_heatmap":
        return [["plot_cn_cell_type_heatmap.py", "--h5ad", cn_aced, cn_bl, "--outdir", final, "--table-dir", tables]]
    if panel == "figure_4bc_whole_tissue":
        return [
            ["analyze_whole_tissue_progression.py", "--analysis", "grade", "--h5ad", cn_aced, cn_bl,
             "--outdir", grade_tables],
            ["plot_whole_tissue_grade.py", "--table-dir", f"{grade_tables}/tables", "--figure-dir", final],
        ]
    if panel == "figure_4e_window_density":
        return [["analyze_axon_window_enrichment.py", "--h5ad", cn_aced, cn_bl,
                 "--outdir", "results/axon_window_enrichment"]]
    if panel == "figure_5b_cn_progression":
        return [
            ["analyze_cn_progression.py", "--group-by", "progression", "--h5ad", cn_aced,
             "--outdir", cn_tables],
            ["plot_cn_progression.py", "--coefficients", f"{cn_tables}/tables/notebook_method_coefficients_aced.csv",
             "--pvalues", f"{cn_tables}/tables/notebook_method_pvalues_aced.csv", "--figdir", final],
        ]
    if panel == "figure_5d_local_redistribution":
        return [["plot_local_redistribution.py", "--master", str(archive_target), "--outdir", final]]
    if panel == "figure_6a_pd1_composition":
        return [
            ["extract_marker_metrics.py", "--input-dir", ppa_dir, "--outdir", "significance_axons_immune_cells",
             "--markers", "PD1"],
            ["plot_pd1_composition.py"],
        ]
    if panel == "figure_6b_whole_tissue_markers":
        return [
            ["analyze_tcell_marker_positivity.py", "--analysis", "grade", "--h5ad", marker_aced, marker_bl,
             "--outdir", marker_grade_tables],
            ["analyze_tcell_marker_positivity.py", "--analysis", "progression", "--h5ad", marker_aced,
             "--outdir", marker_progression_tables],
            ["plot_tcell_marker_grade.py", "--table-dir", f"{marker_grade_tables}/tables", "--figure-dir", final],
            ["plot_tcell_marker_progression.py", "--table-dir", f"{marker_progression_tables}/tables", "--figure-dir", final],
        ]
    if panel == "figure_6cd_axon_marker_fractions":
        return [
            ["extract_marker_metrics.py", "--input-dir", ppa_dir, "--outdir", "significance_axons_immune_cells",
             "--markers", "PD1", "GZMB"],
            ["plot_pd1_progression_barplot.py", "--marker", "PD1"],
            ["plot_pd1_progression_barplot.py", "--marker", "GZMB"],
        ]
    if panel == "figure_6ef_progression_effects":
        return [["plot_progression_marker_effects.py", "--master", str(archive_target), "--outdir", final]]
    raise ValueError(f"Unknown panel: {panel}")


def run_panel(panel: str, data_dir: Path, output_dir: Path) -> None:
    if panel not in PANEL_FILES:
        raise ValueError(f"Unknown panel: {panel}")
    data_dir = data_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if not data_dir.is_dir():
        raise ValueError(f"Data directory does not exist: {data_dir}")

    file_paths = {**H5AD_FILES, **PPA_FILES}
    required_names = [
        "local_progression_statistics.csv" if name == "local_progression_statistics.csv" else file_paths[name]
        for name in PANEL_FILES[panel]
    ]
    inputs = {name: locate_input(data_dir, name) for name in required_names}
    missing = [name for name, path in inputs.items() if path is None]
    if missing:
        raise ValueError("Missing required inputs in data directory:\n  " + "\n  ".join(missing))

    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{panel}_") as temp:
        work = Path(temp)
        environment = os.environ.copy()
        environment.update(
            MPLBACKEND="Agg",
            MPLCONFIGDIR=str(work / "matplotlib"),
            PYTHONHASHSEED="0",
            OMP_NUM_THREADS="1",
            OPENBLAS_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            VECLIB_MAXIMUM_THREADS="1",
        )
        if "local_progression_statistics.csv" in PANEL_FILES[panel]:
            archive_target = work / "manuscript_panels_results/tables/focused_fdr_master_analysis_table.csv"
            archive_target.parent.mkdir(parents=True, exist_ok=True)
            archive_input = inputs["local_progression_statistics.csv"]
            shutil.copyfile(archive_input, archive_target)
            selected_commands = []
        else:
            archive_target = work / "manuscript_panels_results/tables/focused_fdr_master_analysis_table.csv"
            selected_commands = []

        selected_commands.extend(panel_commands(panel, inputs, archive_target))

        for command in selected_commands:
            log = work / f"{Path(command[0]).stem}.log"
            with log.open("w") as handle:
                result = subprocess.run(
                    [sys.executable, str(WORKFLOW / command[0]), *command[1:]],
                    cwd=work,
                    env=environment,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            if result.returncode:
                tail = "\n".join(log.read_text(errors="replace").splitlines()[-30:])
                raise RuntimeError(f"{command[0]} failed:\n{tail}")

        for filename in FIGURE_OUTPUTS[panel]:
            matches = list(work.rglob(filename))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one {filename}, found {len(matches)}")
            source = matches[0]
            for suffix in (".pdf", ".png", ".svg"):
                asset = source.with_suffix(suffix)
                if asset.is_file():
                    shutil.copyfile(asset, output_dir / asset.name)
    print(f"Created {len(FIGURE_OUTPUTS[panel])} figure file(s) in {output_dir}")


def main(panel: str, description: str) -> None:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--data-dir", type=Path, default=Path("data"), help="Flat directory of downloaded Zenodo inputs.")
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    run_panel(panel, args.data_dir, args.output_dir)
