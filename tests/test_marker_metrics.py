"""Small H5 fixtures exercise scientifically important counting/mapping cases."""

from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow"))
from extract_marker_metrics import extract_file, unique_metadata
from patient_ids import explicit_mapping, patient_id_from_region, validate_patient_groups


class MarkerMetricsTests(unittest.TestCase):
    def test_ppa_occurrences_and_zero_denominators(self):
        tag = "CXCR4pos___CD44pos___FYNpos_Axon"
        axon = "CXCR4+ / CD44+ / FYN+ Axon"
        regions = ["/Batch3/001_1_A.ome.tif", "/Batch3/002_1_A.ome.tif"]
        obs = pd.DataFrame({
            "unique_region": [regions[0]] * 4 + [regions[1]] * 3,
            "condition": ["Low-grade"] * 7,
            "progression": ["NonProgressor"] * 4 + ["Progressor"] * 3,
            "cell_type_sync": ["T Cells", "T Cells", "T Cells", axon, "T Cells", "T Cells", axon],
            "cell_type_sync_PD1": ["PD1+ T Cells", "PD1+ T Cells", "PD1- T Cells", axon, "PD1+ T Cells", "PD1- T Cells", axon],
        }, index=[str(i) for i in range(7)])
        data = ad.AnnData(obs=obs)
        # The same positive cell belongs to two stored neighborhoods. Both
        # occurrences must count; zero-neighborhood region remains abundance=0.
        for radius in (20, 50, 100):
            data.uns[f"ppa_result_{radius}"] = obs.iloc[[0, 0, 2]].copy()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.h5ad"
            data.write_h5ad(path)
            composition, abundance = extract_file(path, "PD1", tag)
        self.assertEqual(len(composition), 3)
        self.assertTrue(composition.raw_positive_near.eq(2).all())
        self.assertTrue(composition.raw_marker_split_near_total.eq(3).all())
        self.assertTrue(composition.pct_marker_positive_near.between(66.66, 66.67).all())
        self.assertEqual(len(abundance), 12)
        self.assertTrue(abundance.loc[abundance.unique_region.eq(regions[1]), "raw_n_near"].eq(0).all())

    def test_missing_grade_retains_progression_region(self):
        data = pd.DataFrame({"unique_region": ["r", "r"], "condition": ["Na", "Na"],
                             "progression": ["NonProgressor", "nan"]})
        meta = unique_metadata(data)
        self.assertEqual(meta.loc["r", "progression"], "NonProgressor")
        self.assertEqual(meta.loc["r", "condition"], "Na")

    def test_conflicting_region_metadata_rejected(self):
        data = pd.DataFrame({"unique_region": ["r", "r"], "condition": ["Low-grade", "Low-grade"],
                             "progression": ["NonProgressor", "Progressor"]})
        with self.assertRaises(ValueError):
            unique_metadata(data)

    def test_patient_mapping_merges_samples_not_patient_suffixes(self):
        explicit_mapping.cache_clear()
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(patient_id_from_region("/Batch3/017_1_A_(4)_5.ome.tif"), "017_1")
            self.assertEqual(patient_id_from_region("/Batch3/017_1_B_(4)_5.ome.tif"), "017_1")
            self.assertEqual(patient_id_from_region("/Batch3/017_2_A_(4)_5.ome.tif"), "017_2")

    def test_patient_mapping_conflicting_status_rejected(self):
        explicit_mapping.cache_clear()
        data = pd.DataFrame({"unique_region": ["/Batch3/017_1_A_x", "/Batch3/017_1_B_x"],
                             "progression": ["Progressor", "NonProgressor"]})
        with patch.dict(os.environ, {}, clear=True), self.assertRaises(ValueError):
            validate_patient_groups(data)

    def test_explicit_patient_map_requires_all_regions(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "map.csv"
            path.write_text("unique_region,patient_id\nr1,subject01\n")
            with patch.dict(os.environ, {"AXON_IMMUNE_PATIENT_MAP": str(path)}):
                explicit_mapping.cache_clear()
                self.assertEqual(patient_id_from_region("r1"), "subject01")
                with self.assertRaises(ValueError):
                    patient_id_from_region("r2")
        explicit_mapping.cache_clear()


if __name__ == "__main__":
    unittest.main()
