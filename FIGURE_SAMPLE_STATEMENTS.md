# Figure-specific patient counts and statistical units

Prepared 23 September 2026.

These statements cover all 31 PDFs in the manuscript figure collection. ACED patient IDs follow the current two-token convention (for example, `017_1`), with A/B samples and additional region suffixes nested within patient. For BL, the user confirmed on 23 September 2026 that each numeric `unique_region` identifies a distinct patient. Counts refer to the observations actually available for the displayed endpoint, not the full study enrollment.

The original PD1 donut plot is descriptive: its inferential wording must refer to a separate companion test. Patient-level sensitivity plots pool regions before testing. Archived permutation panels cannot be certified as patient-adjusted. Exact patient counts for reduced radial-distribution subsets remain unverified.

## Statements

### 1. fig_axon_immune_target_window_density_background_vs_axon_aced_bl_combined.pdf

The source dataset includes 18 unique patients contributing 25 regions: 8 ACED patients contributing 15 regions and 10 BL patients contributing one region each. Evaluable comparisons include 18 patients/25 regions for T cells, CD4+ T cells, macrophages, and Tryp+ cells; 15 patients/21 regions for CD20+ cells; 16 patients/23 regions for CD66b+ cells; and 11 patients/14 regions for FOXP3+ T cells. Cell counts and sampled window areas are aggregated within each region to obtain paired axon-associated and background densities. Wilcoxon signed-rank tests use one paired summary per region. Pairing accounts for the two window contexts within a region, but multiple regions from the same patient are treated as independent replicates.

### 2. fig_all_immune_window_density_background_vs_axon_aced_bl_combined.pdf

This analysis includes 18 unique patients contributing 25 regions: 8 ACED patients contributing 15 regions and 10 BL patients contributing one region each. Each region contributes one axon-associated and one background immune-density summary. The paired Wilcoxon signed-rank test uses regions as the unit of analysis; multiple regions from the same patient are treated as independent replicates.

### 3. fig_whole_tissue_grade_composition_separate_panel.pdf

This analysis includes 18 unique patients: 8 from ACED and 10 from BL. The low-grade group comprises 17 regions from 12 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Cell counts are pooled within each grade group to calculate the displayed composition. This panel is descriptive and contains no hypothesis test or p-value.

### 4. fig_whole_tissue_grade_region_immune_fraction_separate_panel.pdf

This analysis includes 18 unique patients: 8 from ACED and 10 from BL. The low-grade group comprises 17 regions from 12 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Cell observations are aggregated within each region to calculate the immune-cell fraction. Groups are compared using a two-sided Mann–Whitney U test on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied. Dependence arising from patients represented in both grade groups is also unaccounted for.

### 5. fig_whole_tissue_grade_immune_density_separate_panel.pdf

This analysis includes 18 unique patients: 8 from ACED and 10 from BL. The low-grade group comprises 17 regions from 12 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Cell observations are aggregated within each region to calculate the immune-cell density using the tissue-area proxy. Groups are compared using a two-sided Mann–Whitney U test on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied. Dependence arising from patients represented in both grade groups is also unaccounted for.

### 6. fig_whole_tissue_tcell_marker_grade_whole_sample_fraction_separate_panel.pdf

This analysis includes 17 unique patients: 7 from ACED and 10 from BL. For each marker, the low-grade group comprises 15 regions from 11 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Marker-positive and classified T-cell counts are pooled across regions within each grade group to calculate the displayed fractions. This panel is descriptive and contains no hypothesis test or p-value.

### 7. fig_whole_tissue_tcell_marker_grade_region_fraction_separate_panel.pdf

This analysis includes 17 unique patients: 7 from ACED and 10 from BL. For each marker, the low-grade group comprises 15 regions from 11 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Cell observations are aggregated within each region to calculate the marker-positive fraction among classified T cells. Groups are compared using two-sided Mann–Whitney U tests on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied. Dependence arising from patients represented in both grade groups is also unaccounted for.

### 8. fig_whole_tissue_tcell_marker_grade_positive_tcell_density_separate_panel.pdf

This analysis includes 17 unique patients: 7 from ACED and 10 from BL. For each marker, the low-grade group comprises 15 regions from 11 patients, and the high-grade group comprises 8 regions from 8 patients. Two ACED patients contribute regions to both grade groups. Cell observations are aggregated within each region to calculate the marker-positive T-cell density using the tissue-area proxy. Groups are compared using two-sided Mann–Whitney U tests on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied. Dependence arising from patients represented in both grade groups is also unaccounted for.

### 9. fig_whole_tissue_tcell_marker_progression_whole_sample_fraction_separate_panel.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients, contributing 15 and 8 regions, respectively. Marker-positive and classified T-cell counts are pooled across regions within each progression group to calculate the displayed fractions. This panel is descriptive and contains no hypothesis test or p-value.

### 10. fig_whole_tissue_tcell_marker_progression_region_fraction_separate_panel.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients, contributing 15 and 8 regions, respectively. Cell observations are aggregated within each region to calculate the marker-positive fraction among classified T cells. Groups are compared using two-sided Mann–Whitney U tests on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied.

### 11. fig_whole_tissue_tcell_marker_progression_positive_tcell_density_separate_panel.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients, contributing 15 and 8 regions, respectively. Cell observations are aggregated within each region to calculate the marker-positive T-cell density using the tissue-area proxy. Groups are compared using two-sided Mann–Whitney U tests on region-level values. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied.

### 12. fig_progression_marker_specific_changes_abundance_separate_panel.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients, contributing 15 and 8 regions, respectively. Nearby AR+, Ki67+, and GZMB+ T-cell counts and tissue totals are summarized within each region for the abundance analysis around NRXN1+/MAOA+ axons. The archived results describe a region-level negative-binomial model with HC0 covariance; HC0 is not a patient-clustered covariance. No patient-level correction is documented, and the original fitting code is unavailable. The Ki67 comparison at 20 µm was not fitted.

### 13. fig_progression_marker_specific_changes_pd1_fraction_separate_panel.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby PD1-positive and PD1-negative T-cell counts are aggregated within each region to calculate one PD1+ fraction per region. The archived results describe a region-level permutation test. Patient-blocked permutation cannot be verified because the original fitting code is unavailable; this panel should not be described as patient-adjusted.

### 14. fig_progression_local_redistribution.pdf

The source cohort comprises 7 unique nonprogressor patients and 5 unique progressor patients. Most displayed cell-type comparisons use all 15 nonprogressor and 8 progressor regions. CD66b+ comparisons use 14 versus 8 regions, and FOXP3+ T-cell comparisons use 10 versus 5 regions. Exact unique-patient counts for the reduced subsets cannot be verified because the archive lacks the contributing region IDs. Spatial observations are summarized within regions to calculate proximal fractions and proximal-to-distal density ratios. The saved analyses are region-level permutation tests; patient-blocked permutation cannot be verified because the original fitting code is unavailable.

### 15. PD1_positive_fraction_barplot.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby PD1-positive and PD1-negative T-cell counts are aggregated within each region to calculate one PD1+ fraction per region. Groups are compared using two-sided Mann–Whitney U tests on region-level fractions. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied.

### 16. GZMB_positive_fraction_barplot.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby GZMB-positive and GZMB-negative T-cell counts are aggregated within each region to calculate one GZMB+ fraction per region. Groups are compared using two-sided Mann–Whitney U tests on region-level fractions. Multiple regions from the same patient are treated as independent region-level observations in the statistical analysis; no patient-level clustering correction is applied.

### 17. PD1_raw_composition_donuts.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby PD1-positive and PD1-negative T-cell counts are pooled across regions within each axon population, radius, and progression group to calculate the displayed fractions. The donuts are descriptive and contain no hypothesis test or p-value. The accompanying region-level comparisons are reported separately.

### 18. PD1_patient_level_composition_donuts.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby PD1-positive and PD1-negative T-cell counts are pooled across regions within each patient to calculate one PD1+ fraction per patient for each axon population and radius. The donuts display the equal-weight mean of these patient fractions. This panel is descriptive and contains no hypothesis test or p-value.

### 19. fig_progression_marker_specific_changes_abundance_patient_level.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients, contributing 15 and 8 regions, respectively. For each marker and radius, nearby marker-positive T-cell counts and tissue totals are pooled across regions within each patient before calculating one tissue-normalized abundance value per patient. Two-sided Mann–Whitney U tests compare the resulting patient-level values. Multiple regions are therefore not treated as independent replicates. This panel uses patient aggregation, not a mixed-effects or cluster-robust model.

### 20. fig_progression_marker_specific_changes_pd1_fraction_patient_level.pdf

This analysis includes 7 unique nonprogressor patients and 5 unique progressor patients. The CXCR4+/CD44+/FYN+ comparisons include 15 nonprogressor and 8 progressor regions at each radius. The NRXN1+/MAOA+ comparisons include 11 versus 7 regions at 20 µm, 12 versus 8 at 50 µm, and 13 versus 8 at 100 µm. All NRXN1+/MAOA+ comparisons include 7 nonprogressor patients; the 20-µm comparison includes 4 progressor patients, whereas the 50- and 100-µm comparisons include 5. Nearby PD1-positive and PD1-negative T-cell counts are pooled across regions within each patient to calculate one PD1+ fraction per patient for each axon population and radius. Two-sided Mann–Whitney U tests compare the patient-level fractions, while plotted effects show differences in patient means. Multiple regions from a patient are combined before testing and are not treated as independent replicates.

### 21. aced_progression_paper_ready.pdf

This figure includes 8 unique ACED patients: 3 nonprogressors contributing 7 regions and 5 progressors contributing 8 regions. Regions sharing a biopsy/time-point identifier are pooled, producing 5 nonprogressor and 7 progressor samples for the CN analysis. Cell and axon observations are aggregated within these samples to calculate CN-specific fractions. The OLS models use the 12 samples as observations; separate samples from the same patient are treated as independent, without a patient random effect or patient-clustered covariance.

### 22. fig_cn_cell_type_heatmap_aced.pdf

This figure includes 8 unique ACED patients: 3 nonprogressors contributing 7 regions and 5 progressors contributing 8 regions. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 23. fig_cn_cell_type_heatmap_bl.pdf

This figure includes 10 unique BL patients, each contributing one region. Seven patients have low-grade and three have high-grade lesions. Progression status is not available for this cohort. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 24. fig_cn_cell_type_heatmap_combined.pdf

This figure includes 18 unique patients contributing 25 regions: 8 ACED patients contributing 15 regions and 10 BL patients contributing one region each. The ACED subset comprises 3 nonprogressor and 5 progressor patients; progression status is not available for BL. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 25. fig_cn_cell_type_heatmap_cn_rows_aced.pdf

This figure includes 8 unique ACED patients: 3 nonprogressors contributing 7 regions and 5 progressors contributing 8 regions. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 26. fig_cn_cell_type_heatmap_cn_rows_bl.pdf

This figure includes 10 unique BL patients, each contributing one region. Seven patients have low-grade and three have high-grade lesions. Progression status is not available for this cohort. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 27. fig_cn_cell_type_heatmap_cn_rows_combined.pdf

This figure includes 18 unique patients contributing 25 regions: 8 ACED patients contributing 15 regions and 10 BL patients contributing one region each. The ACED subset comprises 3 nonprogressor and 5 progressor patients; progression status is not available for BL. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 28. fig_cn_exp_heatmap_function_aced.pdf

This figure includes 8 unique ACED patients: 3 nonprogressors contributing 7 regions and 5 progressors contributing 8 regions. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 29. fig_cn_exp_heatmap_function_bl.pdf

This figure includes 10 unique BL patients, each contributing one region. Seven patients have low-grade and three have high-grade lesions. Progression status is not available for this cohort. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 30. fig_cn_exp_heatmap_function_combined.pdf

This figure includes 18 unique patients contributing 25 regions: 8 ACED patients contributing 15 regions and 10 BL patients contributing one region each. The ACED subset comprises 3 nonprogressor and 5 progressor patients; progression status is not available for BL. Cell and axon counts are pooled across regions within each cellular neighborhood to calculate the displayed enrichment. This heatmap is descriptive and contains no hypothesis test or p-value. Any hierarchical clustering of heatmap rows or columns is a visualization step, not patient-level statistical adjustment.

### 31. fig_cn_exp_heatmap_color_scale_legend.pdf

This figure is a color-scale legend and contains no patient observations, hypothesis tests, or p-values.

## Window-density comparison counts

| target | n_patients | n_regions |
| --- | --- | --- |
| CD20+ | 15 | 21 |
| CD4+ T Cells | 18 | 25 |
| CD66b+ | 16 | 23 |
| FOXP3+ T Cells | 11 | 14 |
| Macrophages | 18 | 25 |
| T Cells | 18 | 25 |
| Tryp+ | 18 | 25 |

These counts describe paired observations: each contributing region supplies both an axon-associated and a background density. They are not separate sets of patients for the two window contexts.

## Source verification

Counts were calculated from the saved region-level marker-composition, marker-abundance, whole-tissue, and window-summary tables. CN cohort coverage was checked against H5 observation metadata; CN sample counts use the actual A/B sample parser. Statistical units follow the generating code and the prior figure audit. Reduced radial subset region counts come from the archived master table; unique-patient membership is not present there.

`manifests/figure_sample_counts.csv` provides the comparison-specific count register. For descriptive figures, there is no statistical replicate count. Blank patient counts in the reduced radial subsets mean unverified, not zero.
