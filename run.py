#!/usr/bin/env python3
"""Reproduce the manuscript figure collection from H5AD inputs and archived fits."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
WORKFLOW = ROOT / "workflow"
STAGES = ("whole_tissue", "markers", "patients", "windows", "cn", "archived")
CN_ACED = "data/adata_aced_immune_all_CN_k50_n10_with_epithelia_filtered.h5ad"
CN_BL = "data/adata_BLimmune_all_CN_k50_n10_with_epithelia.h5ad"
MARKER_ACED = "data/adata_aced_immune_updated_tcell_marker_celltypes.h5ad"
MARKER_BL = "data/adata_bl_immune_updated_tcell_marker_celltypes.h5ad"
TAGS = ("CXCR4pos___CD44pos___FYNpos_Axon", "NRXN1pos___MAOApos_Axon")
PPA_FILES = [f"PPA_axon_subtypes_marker_split_T_cells/aced_immune_{m}_{a}_ppa.h5ad"
             for m in ("AR", "Ki67", "GZMB", "PD1") for a in TAGS]
STAGE_INPUTS = {
    "whole_tissue": [CN_ACED, CN_BL, MARKER_ACED, MARKER_BL],
    "markers": PPA_FILES, "patients": PPA_FILES,
    "windows": [CN_ACED, CN_BL], "cn": [CN_ACED, CN_BL], "archived": [],
}


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / "inputs", help="Directory containing data/ and PPA_axon_subtypes_marker_split_T_cells/.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    parser.add_argument("--stages", nargs="+", choices=("all",) + STAGES, default=["all"])
    parser.add_argument("--patient-map", type=Path, help="Optional CSV with unique_region,patient_id; otherwise use first two ACED filename tokens.")
    parser.add_argument("--resume", action="store_true", help="Reuse completed stages only when inputs, options, and code match.")
    parser.add_argument("--verify-inputs", action="store_true", help="Verify complete H5AD SHA256 values against the input manifest.")
    parser.add_argument("--check-inputs", action="store_true", help="Validate input layout/metadata and exit.")
    parser.add_argument("--plan", action="store_true", help="Print stages and commands without reading input data or running analyses.")
    return parser.parse_args()


def selected_stages(requested: list[str]) -> list[str]:
    wanted = set(STAGES if "all" in requested else requested)
    if "patients" in wanted:
        wanted.add("markers")
    return [stage for stage in STAGES if stage in wanted]


def stage_commands(data_root: Path) -> dict[str, list[list[str]]]:
    source = lambda relative: str(data_root / relative)
    final = "manuscript final figures/figures"
    tables = "manuscript final figures/tables"
    immune_dir = "manuscript_panels_results/whole_tissue_grade/combined"
    grade_dir = "manuscript_panels_results/whole_tissue_tcell_marker_grade_updated/combined"
    progression_dir = "manuscript_panels_results/whole_tissue_tcell_marker_progression_updated"
    cn_dir = "results/cn_neighborhood_enrichment_notebook_method_progression_aced"
    return {
        "whole_tissue": [
            ["analyze_whole_tissue_progression.py", "--analysis", "grade", "--h5ad", source(CN_ACED), source(CN_BL), "--outdir", immune_dir],
            ["analyze_tcell_marker_positivity.py", "--analysis", "grade", "--h5ad", source(MARKER_ACED), source(MARKER_BL), "--outdir", grade_dir],
            ["analyze_tcell_marker_positivity.py", "--analysis", "progression", "--h5ad", source(MARKER_ACED), "--outdir", progression_dir],
            ["replot_whole_tissue_grade_separate_panels.py", "--table-dir", immune_dir + "/tables", "--figure-dir", final],
            ["replot_tcell_marker_grade_separate_panels.py", "--table-dir", grade_dir + "/tables", "--figure-dir", final],
            ["replot_tcell_marker_progression_separate_panels.py", "--table-dir", progression_dir + "/tables", "--figure-dir", final],
        ],
        "markers": [
            ["extract_marker_metrics.py", "--input-dir", source("PPA_axon_subtypes_marker_split_T_cells")],
            ["replot_pd1_progression_barplot_paper_ready.py", "--marker", "PD1"],
            ["replot_pd1_progression_barplot_paper_ready.py", "--marker", "GZMB"],
            ["replot_pd1_progression_donuts_paper_ready.py"],
        ],
        "patients": [["replot_pd1_patient_level_donuts.py"], ["replot_progression_marker_specific_changes_patient_level.py"]],
        "windows": [["analyze_axon_immune_spatial_enrichment_combined_plot.py", "--h5ad", source(CN_ACED), source(CN_BL),
                     "--outdir", "results/axon_immune_spatial_enrichment_aced_bl_combined_plot"]],
        "cn": [
            ["analyze_cn_neighborhood_enrichment_notebook_method.py", "--group-by", "progression", "--h5ad", source(CN_ACED), "--outdir", cn_dir],
            ["replot_cn_neighborhood_progression_aced_paper_ready.py", "--coefficients", cn_dir + "/tables/notebook_method_coefficients_aced.csv",
             "--pvalues", cn_dir + "/tables/notebook_method_pvalues_aced.csv", "--figdir", final],
            ["plot_cn_cell_type_heatmap_paper_ready.py", "--h5ad", source(CN_ACED), source(CN_BL), "--outdir", final, "--table-dir", tables],
            ["plot_cn_cell_type_heatmap_paper_ready.py", "--h5ad", source(CN_ACED), source(CN_BL), "--outdir", final, "--table-dir", tables,
             "--orientation", "cn_rows", "--prefix", "fig_cn_cell_type_heatmap_cn_rows"],
            ["plot_cn_exp_heatmap_function_paper_ready.py", "--h5ad", source(CN_ACED), source(CN_BL), "--outdir", final, "--table-dir", tables],
            ["plot_cn_heatmap_color_scale_legend.py", "--outdir", final],
        ],
        "archived": [
            ["replot_progression_marker_specific_changes_separate_panels.py", "--outdir", final],
            ["replot_progression_local_redistribution_transposed.py", "--outdir", final],
        ],
    }


def check_inputs(data_root: Path, stages: list[str], verify: bool) -> list[dict]:
    required = sorted({name for stage in stages for name in STAGE_INPUTS[stage]})
    missing = [name for name in required if not (data_root / name).is_file()]
    if missing:
        raise ValueError("Missing input files under " + str(data_root) + ":\n  " + "\n  ".join(missing))
    if not required:
        return []
    import h5py
    expected = {}
    manifest = ROOT / "manifests/inputs.csv"
    if manifest.exists():
        with manifest.open(newline="") as handle:
            expected = {row["relative_path"]: row for row in csv.DictReader(handle)}
    records = []
    for name in required:
        path = data_root / name
        with h5py.File(path, "r") as handle:
            required_columns = {"unique_region", "cell_type_sync", "condition"}
            if "aced" in name.lower():
                required_columns.add("progression")
            if name.startswith("data/"):
                required_columns |= {"x", "y"}
            if name in (CN_ACED, CN_BL):
                required_columns.add("CN_k50_n10")
            if "obs" not in handle or not required_columns.issubset(handle["obs"].keys()):
                raise ValueError(f"Required observation metadata missing from {name}")
            if name in PPA_FILES:
                marker = path.name.split("_")[2]
                for radius in (20, 50, 100):
                    key = f"uns/ppa_result_{radius}"
                    if key not in handle or f"cell_type_sync_{marker}" not in handle[key]:
                        raise ValueError(f"Required PPA result missing from {name}: {key}")
        record = dict(relative_path=name, bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
        if verify:
            print(f"Checking SHA256: {name}", flush=True)
            record["sha256"] = sha256(path)
            if name not in expected or record["sha256"] != expected[name]["sha256"]:
                raise ValueError(f"Input differs from the reference H5AD: {name}")
        records.append(record)
    return records


def collect_figures(stage: str, work: Path, destination: Path, manifest: list[dict]) -> list[dict]:
    destination.mkdir(parents=True, exist_ok=True)
    records = []
    for row in manifest:
        if row["stage"] != stage:
            continue
        matches = list(work.rglob(row["filename"]))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one generated {row['filename']}, found {len(matches)}")
        source = matches[0]
        for suffix in (".pdf", ".png", ".svg"):
            other = source.with_suffix(suffix)
            if other.exists():
                shutil.copyfile(other, destination / other.name)
        records.append(dict(row, generated_from=str(source.relative_to(work)), sha256=sha256(destination / source.name)))
    return records


def package_versions() -> dict[str, str]:
    result = {"python": sys.version.split()[0]}
    for line in (ROOT / "requirements.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name = line.split("==")[0]
            result[name] = importlib.metadata.version(name)
    return result


def main() -> None:
    args = arguments()
    stages = selected_stages(args.stages)
    data_root, output = args.data_root.resolve(), args.output_dir.resolve()
    commands = stage_commands(data_root)
    if args.plan:
        for stage in stages:
            print(stage + ":")
            for command in commands[stage]:
                print("  " + repr(command))
        return
    if output == ROOT or output == data_root or data_root in output.parents and output.name in ("data", "PPA_axon_subtypes_marker_split_T_cells"):
        raise ValueError("Use a dedicated output directory, separate from the repository and input files.")
    inputs = check_inputs(data_root, stages, args.verify_inputs)
    if args.patient_map and not args.patient_map.is_file():
        raise ValueError(f"Patient map does not exist: {args.patient_map}")
    if args.check_inputs:
        print(f"Validated {len(inputs)} H5AD files for stages: {', '.join(stages)}")
        return
    versions = package_versions()
    output.mkdir(parents=True, exist_ok=True)
    work = output / "analysis"
    logs = output / "logs"
    work.mkdir(exist_ok=True)
    logs.mkdir(exist_ok=True)
    env = os.environ.copy()
    env.update(MPLBACKEND="Agg", MPLCONFIGDIR=str(logs / "matplotlib"), PYTHONHASHSEED="0", PYTHONUNBUFFERED="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    if args.patient_map:
        env["AXON_IMMUNE_PATIENT_MAP"] = str(args.patient_map.resolve())
    else:
        env.pop("AXON_IMMUNE_PATIENT_MAP", None)
    with (ROOT / "manifests/figures.csv").open(newline="") as handle:
        figures = list(csv.DictReader(handle))
    code_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in sorted(WORKFLOW.glob("*.py"))}
    code_hashes["run.py"] = sha256(Path(__file__))
    base = dict(data_root=str(data_root), inputs=inputs, versions=versions, code=code_hashes,
                patient_map_sha256=sha256(args.patient_map) if args.patient_map else None,
                patient_rule="explicit map" if args.patient_map else "first two ACED filename tokens",
                archived_sha256=sha256(ROOT / "archived/local_progression_statistics.csv"))
    fingerprint = hashlib.sha256(json.dumps(base, sort_keys=True).encode()).hexdigest()
    provenance = dict(base, started_at=utc_now(), requested_stages=stages, stages=[])
    all_figures = []
    for stage in stages:
        status_file = logs / f"{stage}.completed.json"
        if args.resume and status_file.exists():
            previous = json.loads(status_file.read_text())
            expected_figures = [row["filename"] for row in figures if row["stage"] == stage]
            if previous.get("fingerprint") == fingerprint and all((output / "figures" / name).exists() for name in expected_figures):
                print(f"[{stage}] reusing verified completed stage", flush=True)
                all_figures.extend(collect_figures(stage, work, output / "figures", figures))
                provenance["stages"].append(dict(stage=stage, status="reused"))
                continue
        if status_file.exists():
            status_file.unlink()
        if stage == "archived":
            archive_target = work / "manuscript_panels_results/tables/focused_fdr_master_analysis_table.csv"
            archive_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / "archived/local_progression_statistics.csv", archive_target)
        start = time.monotonic()
        for number, command in enumerate(commands[stage], start=1):
            log = logs / f"{stage}_{number:02d}_{Path(command[0]).stem}.log"
            print(f"[{stage} {number}/{len(commands[stage])}] {command[0]}", flush=True)
            with log.open("w") as handle:
                result = subprocess.run([sys.executable, str(WORKFLOW / command[0]), *command[1:]], cwd=work, env=env,
                                        stdout=handle, stderr=subprocess.STDOUT, check=False)
            if result.returncode:
                tail = "\n".join(log.read_text(errors="replace").splitlines()[-30:])
                raise RuntimeError(f"{command[0]} failed; see {log}\n{tail}")
        generated = collect_figures(stage, work, output / "figures", figures)
        all_figures.extend(generated)
        status = dict(stage=stage, status="completed", seconds=round(time.monotonic() - start, 2),
                      fingerprint=fingerprint, completed_at=utc_now(), commands=commands[stage])
        status_file.write_text(json.dumps(status, indent=2) + "\n")
        provenance["stages"].append(status)
        (output / "run_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(f"[{stage}] {len(generated)} figures complete ({status['seconds']} seconds)", flush=True)
    with (output / "generated_figures.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_figures[0]))
        writer.writeheader()
        writer.writerows(all_figures)
    provenance["finished_at"] = utc_now()
    (output / "run_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Complete: {len(all_figures)} figures in {output / 'figures'}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, FileNotFoundError, importlib.metadata.PackageNotFoundError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
