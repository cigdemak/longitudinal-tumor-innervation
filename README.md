# Manuscript code

Code accompanying Ak, Ait-Ahmad et al., “Longitudinal Analysis of Tumor Innervation Identifies Dynamic Neuroimmune Niches in Prostate Cancer,” submitted to *Nature Cancer*.

Annotated data tables are available from [Zenodo](https://zenodo.org/records/22968980).
Use [zenodo_upload_inputs.txt](zenodo_upload_inputs.txt) to collect the files needed to reproduce the listed figure panels.

## Run a figure script

Use Python 3.12. Create the environment once, then run any script from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python figure_code/figure_4a_cn_heatmap.py --data-dir /path/to/zenodo-data
```

For the existing sibling-folder layout, run with `--data-dir ..` from the repository root.

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

See [zenodo_upload_inputs.txt](zenodo_upload_inputs.txt) for the required files.
