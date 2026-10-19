# Manuscript code

Code accompanying Ak, Ait-Ahmad et al., “Longitudinal Analysis of Tumor Innervation Identifies Dynamic Neuroimmune Niches in Prostate Cancer,” submitted to *Nature Cancer*.

The H5AD inputs needed to reproduce the figure panels are listed in
[zenodo_upload_inputs.txt](zenodo_upload_inputs.txt). The analysis scripts read
the cell-level observations and PPA neighborhood tables stored in those files.
The Zenodo-facing filenames use the manuscript cohort labels AS (tumor) and RP
(benign/reference prostate); the runner also recognizes the original analysis
filenames in the local project folders.

## Run a figure script

Use Python 3.12. Create the environment once, then run all panels listed below
from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python figure_code/reproduce_all.py --data-dir /path/to/h5ad-files
```

Figures are written to `figures/`. To reproduce one panel, run its
script directly, for example:

```bash
python figure_code/figure_5d_local_redistribution.py \
  --data-dir /path/to/h5ad-files \
  --output-dir figures
```

## Figure scripts

- `figure_4a_cn_heatmap.py` (Figure 4a)
- `figure_4bc_whole_tissue.py` (Figure 4b–c)
- `figure_4e_window_density.py` (Figure 4e)
- `figure_5b_cn_progression.py` (Figure 5b)
- `figure_5d_local_redistribution.py` (Figure 5d)
- `figure_6a_pd1_composition.py` (Figure 6a)
- `figure_6b_whole_tissue_markers.py` (Figure 6b)
- `figure_6cd_axon_marker_fractions.py` (Figure 6c–d)
- `figure_6ef_progression_effects.py` (Figure 6e–f)

See [zenodo_upload_inputs.txt](zenodo_upload_inputs.txt) for the complete H5AD
input list. Install dependencies from `requirements.txt` before running.
