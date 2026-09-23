# Axon–immune manuscript figures

Reproduce the 31 figures in the manuscript figure collection from annotated
H5AD files. **The H5 workflow generates 28 figures (including the color-scale
legend); three original effect panels use included archived statistics.**

## Run

Use Python 3.12:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run.py --data-root /path/to/input-folder
```

The input folder must contain:

```text
input-folder/
├── data/                                       # four cohort H5AD files
└── PPA_axon_subtypes_marker_split_T_cells/      # eight ACED marker/axon H5AD files
```

Exact filenames, sizes, and SHA256 checksums are in
[manifests/inputs.csv](manifests/inputs.csv). H5 files are supplied separately,
not included in this repository. Inputs already contain cell annotations,
CN labels, and `uns/ppa_result_20`, `_50`, and `_100` neighborhoods; this
repository starts from those annotations and stored neighborhoods.

Outputs:

- `outputs/figures/`: PDF/PNG figures and SVG where supported.
- `outputs/analysis/`: regenerated analysis tables and intermediate exports.
- `outputs/generated_figures.csv`: figure provenance and statistical units.
- `outputs/logs/` and `outputs/run_provenance.json`: commands and environment.

## Options

```bash
python run.py --data-root /path/to/input-folder --check-inputs --verify-inputs
python run.py --data-root /path/to/input-folder --resume
python run.py --data-root /path/to/input-folder --stages markers patients
python run.py --stages archived
python -m unittest discover -s tests -v
python verify_figures.py
```

Stages: `whole_tissue`, `markers`, `patients`, `windows`, `cn`, `archived`.
Use `--output-dir` to choose the output folder; `--plan` prints the commands.
The `patients` stage automatically includes `markers`.
The full run reproduced all 31 reference PNGs exactly; see [VALIDATION.md](VALIDATION.md).

## Reproduction scope

[manifests/figures.csv](manifests/figures.csv) lists each figure, its stage,
statistical unit, and p/q annotation. Existing statistical choices are
preserved, including region/sample-level tests without patient clustering.

Patient sensitivity figures pool ACED regions using the first two filename
tokens as the patient ID (`016_1`). To supply a different verified mapping,
use `--patient-map mapping.csv` with columns `unique_region,patient_id`.
Conflicting progression labels within a mapped patient cause an error.
The original CN model retains its separate A/B sample IDs.
For BL, each numeric `unique_region` identifies one patient.
[Figure-specific statements](FIGURE_SAMPLE_STATEMENTS.md) report patient and
region counts and how repeated observations are handled for all 31 figures;
[comparison-specific counts](manifests/figure_sample_counts.csv) are also included.

The local abundance, PD1 progression-effect, and radial-redistribution panels
use [archived/local_progression_statistics.csv](archived/local_progression_statistics.csv).
Their original fitting code is unavailable; these three model fits are **not
recomputed from H5**. The archive retains the full correction families.
The other local marker counts are reconstructed from stored PPA rows,
including overlapping cell–axon-neighborhood occurrences.
