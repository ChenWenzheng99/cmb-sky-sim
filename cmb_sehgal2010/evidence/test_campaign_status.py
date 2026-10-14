"""Bounded-I/O status helper controls; no native map creation or pixel reads."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("campaign_status_test", BASE / "campaign_status.py")
status = importlib.util.module_from_spec(spec)
spec.loader.exec_module(status)
PRODUCER_PATH = BASE.parent / "production" / "lensed_CMB_generation.py"


class CampaignStatusTests(unittest.TestCase):
    def test_capture_and_stat_changes_without_reading_pixel_contents(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for rid in [0, 1]:
                paths = status.paths_for(root, rid)
                paths["directory"].mkdir()
                paths["old"].write_bytes(b"not even a FITS; stat only")
                paths["kappa"].write_bytes(b"stat-only kappa")
            inventory = status.capture_inventory(root, [0, 1])
            self.assertEqual(inventory["old_cmb_count"], 2)
            self.assertEqual(inventory["kappa_count"], 2)
            self.assertTrue(inventory["source_filename_counts_exact"])
            self.assertFalse(inventory["large_source_files_hashed"])
            old = inventory["realizations"][0]["old_cmb"]
            status.verify_stat_identity(status.stat_record(Path(old["path"])), old, "old")
            Path(old["path"]).write_bytes(b"changed")
            with self.assertRaises(ValueError):
                status.verify_stat_identity(status.stat_record(Path(old["path"])), old, "old")
            status.paths_for(root, 0)["directory"].joinpath("kappa_00001.fits").write_bytes(b"extra")
            extra = status.capture_inventory(root, [0, 1])
            self.assertEqual(extra["kappa_count"], 3)
            self.assertFalse(extra["source_filename_counts_exact"])

    def test_expected_seed_independent_of_producer_implementation(self):
        producer, _ = status.load_producer(PRODUCER_PATH)
        for rid in [0, 1, 9, 499]:
            self.assertEqual(status.independent_seed(20261011, rid), producer.realization_seed(20261011, rid))
        self.assertEqual(len({status.independent_seed(20261011, rid) for rid in range(500)}), 500)

    def test_finite_record_rejects_missing_unseen_and_mismatched_serialization(self):
        field = {"shape": [status.NPIX], "dtype": "float64", "count": status.NPIX,
                 "finite_count": status.NPIX, "unseen_count": 0,
                 "min": -1., "max": 2., "mean": .2, "population_std": 1.}
        metadata = {"finite_statistics": {"source_kappa": field,
                   "generated_TQU": [dict(field) for _ in range(3)],
                   "output_serialized_TQU": [dict(field) for _ in range(3)]}}
        status.verify_finite_statistics(metadata)
        metadata["finite_statistics"]["output_serialized_TQU"][1]["unseen_count"] = 1
        with self.assertRaises(ValueError):
            status.verify_finite_statistics(metadata)

    def test_native_header_contract_and_checksum_presence_only(self):
        producer, _ = status.load_producer(PRODUCER_PATH)
        metadata = {"coordinate_provenance": {"output_coordinate": "G"},
                    "theory_archive": {"sha256": "a" * 64}, "source_kappa": {"sha256": "b" * 64}}
        header = {"NSIDE": 2048, "ORDERING": "RING", "COORDSYS": "G", "POLCCONV": "COSMO",
                  "TFIELDS": 3, "INDXSCHM": "IMPLICIT", "FIRSTPIX": 0, "LASTPIX": status.NPIX - 1,
                  "REALID": 0, "RNGSEED": 15, "THSHA": "a" * 64, "KPSHA": "b" * 64,
                  "UNLLMAX": 5999, "DEFLMAX": 6143, "KAPITER": 3, "REMEPS": 1e-6,
                  "NAXIS1": 3 * 1024 * 8, "NAXIS2": status.NPIX // 1024,
                  **{f"TTYPE{i}": name for i, name in enumerate(["T", "Q", "U"], 1)},
                  **{f"TUNIT{i}": "uK_CMB" for i in [1, 2, 3]},
                  **{f"TFORM{i}": "1024D" for i in [1, 2, 3]},
                  "checksum_cards": [{"CHECKSUM": "a" * 16, "DATASUM": "0"}]}
        status.verify_header(header, metadata, 0, 15, producer)
        for key, bad in [("NSIDE", 512), ("KAPITER", 1), ("NAXIS2", 1), ("RNGSEED", 16), ("COORDSYS", "C")]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                status.verify_header({**header, key: bad}, metadata, 0, 15, producer)
        with self.assertRaises(ValueError):
            status.verify_header({**header, "checksum_cards": [{"CHECKSUM": None, "DATASUM": None}]},
                                 metadata, 0, 15, producer)


if __name__ == "__main__":
    unittest.main(verbosity=2)
