"""Small controls of the producer contract, without native-map generation."""
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import healpy as hp
import numpy as np

MODULE_PATH = Path(__file__).with_name("lensed_CMB_generation.py")
spec = importlib.util.spec_from_file_location("sehgal_production", MODULE_PATH)
producer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(producer)


def synthetic_theory(directory, change=None):
    """Minimal on-disk fixture obeying the actual fixed-cutoff input contract."""
    ell = np.arange(6501, dtype=np.int64)
    primary = np.zeros((6501, 4), dtype=np.float64)
    primary[2:, 0] = 2 / (ell[2:] * (ell[2:] + 1))
    primary[2:, 1] = .4 / (ell[2:] * (ell[2:] + 1))
    primary[2:, 3] = .1 / (ell[2:] * (ell[2:] + 1))
    cut = primary.copy()
    cut[6000:] = 0
    phi = np.zeros(6501, dtype=np.float64)
    phi[2:] = 1e-10 / (ell[2:] * (ell[2:] + 1)) ** 2
    arrays = {"ell": ell, "unlensed_raw_cl": primary, "lensed_raw_cl": primary.copy(),
              "producer_unlensed_raw_cl": cut, "phi_raw_cl": phi,
              "phi_PP": phi * (ell * (ell + 1)) ** 2 / (2 * np.pi)}
    if change:
        change(arrays)
    archive = directory / "sehgal2010_theory.npz"
    np.savez_compressed(archive, **arrays)
    manifest = {"schema": producer.THEORY_SCHEMA,
                "status": "COMPLETE_THEORY_WITH_DECLARED_ASSUMPTIONS",
                "producer_cutoffs": {"primary_lmax": 5999, "kappa_deflection_lmax": 6143, "nside": 2048},
                "array_contract": {key: {"shape": list(array.shape), "dtype": str(array.dtype),
                                   "units": "uK_CMB^2" if array.ndim == 2 else "dimensionless",
                                   **({"columns": producer.COLUMNS} if array.ndim == 2 else {})}
                                   for key, array in arrays.items()},
                "artifact_manifest": {archive.name: producer.file_record(archive)}}
    path = directory / "theory_manifest.json"
    path.write_text(json.dumps(manifest))
    return archive, path


class ProductionControlTests(unittest.TestCase):
    def test_original_snapshot_immutable_hash(self):
        records = producer.production_source_records()
        self.assertEqual(records["original_snapshot"]["sha256"], producer.ORIGINAL_SHA256)
        self.assertEqual(records["original_snapshot"]["bytes"], 3000)
        self.assertFalse(Path(records["original_snapshot"]["path"]).stat().st_mode & 0o222)

    def test_seed_order_rank_independence_and_bounds(self):
        expected = {i: producer.realization_seed(20261011, i) for i in range(10)}
        self.assertEqual(expected, {i: producer.realization_seed(20261011, i) for i in reversed(range(10))})
        self.assertEqual(expected, {i: producer.realization_seed(20261011, i) for rank in range(5) for i in range(rank, 10, 5)})
        self.assertEqual(len(set(expected.values())), 10)
        self.assertTrue(all(0 <= seed < 2 ** 32 for seed in expected.values()))
        self.assertNotEqual(expected[0], producer.realization_seed(20261012, 0))
        for value in [-1, True, 2 ** 64]:
            with self.assertRaises(ValueError):
                producer.realization_seed(value, 0)

    def test_explicit_synalm_order_matches_legacy_and_restores_rng(self):
        ell = np.arange(17, dtype=float)
        raw = np.zeros((17, 4))
        raw[2:, 0] = 2 / (ell[2:] * (ell[2:] + 1))
        raw[2:, 1] = .4 / (ell[2:] * (ell[2:] + 1))
        raw[2:, 3] = .1 / (ell[2:] * (ell[2:] + 1))
        np.random.seed(123)
        expected_rng_tail = np.random.random(5)
        np.random.seed(123)
        actual = producer.draw_unlensed_alms(raw, 456, lmax=16)
        np.testing.assert_array_equal(np.random.random(5), expected_rng_tail)
        np.random.seed(456)
        legacy = hp.synalm([raw[:, 0], raw[:, 3], raw[:, 1], raw[:, 2]], lmax=16, mmax=16, new=False)
        np.testing.assert_array_equal(actual, legacy)
        np.testing.assert_array_equal(actual[2], np.zeros_like(actual[2]))
        np.testing.assert_array_equal(actual, producer.draw_unlensed_alms(raw, 456, lmax=16))

    def test_harmonic_phi_and_deflection_sign_factors(self):
        phi, gradient = producer.harmonic_transfers(16)
        ell = np.arange(17)
        self.assertEqual(phi[0], 0)
        self.assertEqual(gradient[0], 0)
        self.assertEqual(phi[1], 1)
        self.assertAlmostEqual(gradient[1], np.sqrt(2))
        np.testing.assert_allclose(phi[1:], 2 / (ell[1:] * (ell[1:] + 1)), rtol=0, atol=0)
        kappa = np.arange(hp.Alm.getsize(16), dtype=float).astype(complex)
        kappa *= 1 + .3j
        kphi = hp.almxfl(kappa, phi)
        kd = hp.almxfl(kappa, gradient)
        cl = hp.alm2cl(kappa)
        np.testing.assert_allclose(hp.alm2cl(kphi), cl * phi ** 2, rtol=1e-15)
        np.testing.assert_allclose(hp.alm2cl(kd), cl * gradient ** 2, rtol=1e-15)
        np.testing.assert_allclose(hp.alm2cl(kappa, kphi), cl * phi, rtol=1e-15)
        np.testing.assert_allclose(hp.alm2cl(kphi, kd), cl * phi * gradient, rtol=1e-15)

    def test_coordinate_header_priority_and_no_rotation(self):
        owner = producer.coordinate_provenance({"COORDSYS": None}, "G")
        self.assertEqual(owner["output_coordinate"], "G")
        self.assertFalse(owner["source_frame_available_in_header"])
        self.assertIn("ICRS", owner["caveat"])
        header = producer.coordinate_provenance({"COORDSYS": "C"}, "G")
        self.assertEqual(header["output_coordinate"], "C")
        self.assertTrue(header["source_frame_available_in_header"])
        self.assertFalse(header["coordinate_rotation_performed"])
        with self.assertRaises(ValueError):
            producer.coordinate_provenance({"COORDSYS": "unknown"})

    def test_full_statistics_finite_unseen_and_byte_order(self):
        values = np.arange(64, dtype=np.float64)
        self.assertEqual(producer.full_statistics(values), producer.full_statistics(values.astype(">f8")))
        self.assertEqual(producer.full_statistics(values)["finite_count"], 64)
        for bad in [np.nan, np.inf, hp.UNSEEN]:
            corrupted = values.copy()
            corrupted[31] = bad
            with self.assertRaises(ValueError):
                producer.full_statistics(corrupted)
        with self.assertRaises(ValueError):
            producer.full_statistics(values, expected_shape=(63,))

    def test_theory_loader_valid_and_rejects_hash_column_mismatch(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            archive, manifest_path = synthetic_theory(directory)
            result = producer.load_theory(archive, manifest_path)
            self.assertEqual(result["arrays"]["unlensed_raw_cl"].shape, (6501, 4))
            self.assertEqual(result["archive"]["sha256"], producer.sha256(archive))
            manifest = json.loads(manifest_path.read_text())
            manifest["array_contract"]["unlensed_raw_cl"]["columns"] = ["TT", "TE", "EE", "BB"]
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "TT/EE/BB/TE"):
                producer.load_theory(archive, manifest_path)
            archive.write_bytes(archive.read_bytes() + b"changed")
            with self.assertRaisesRegex(ValueError, "hash/size"):
                producer.load_theory(archive, manifest_path)

    def test_theory_rejects_unlensed_bb_cutoff_phi_and_ell_errors(self):
        changes = [lambda a: a["unlensed_raw_cl"].__setitem__((5, 2), .01),
                   lambda a: a["producer_unlensed_raw_cl"].__setitem__((6001, 0), .01),
                   lambda a: a["phi_PP"].__setitem__(100, .5),
                   lambda a: a["ell"].__setitem__(100, 99),
                   lambda a: a["unlensed_raw_cl"].__setitem__((5, 3), 10)]
        for change in changes:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temp:
                archive, manifest = synthetic_theory(Path(temp), change)
                with self.assertRaises(ValueError):
                    producer.load_theory(archive, manifest)

    def test_ids_and_immutable_output_names(self):
        self.assertEqual(producer.parse_ids(ids="0,2,499"), [0, 2, 499])
        self.assertEqual(producer.parse_ids(start_id=3, count=2), [3, 4])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "ids.txt"
            path.write_text("0, 2 # note\n499\n")
            self.assertEqual(producer.parse_ids(ids_file=path), [0, 2, 499])
        for options in [{"ids": "0,0"}, {"ids": ""}, {"ids": "0", "start_id": 0, "count": 1},
                        {"start_id": 0}, {"ids": "-1"}]:
            with self.assertRaises(ValueError):
                producer.parse_ids(**options)
        paths = producer.output_paths("/tmp/example", 5)
        self.assertEqual(paths["fits"].name, "lensedcmb_new_00005.fits")
        self.assertEqual(paths["old_cmb"].name, "lensedcmb_00005.fits")
        self.assertNotEqual(paths["fits"], paths["old_cmb"])

    def test_atomic_no_clobber_backup_and_lock(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp) / "00000"
            directory.mkdir()
            paths = producer.output_paths(temp, 0)
            paths["old_cmb"].write_bytes(b"old map must remain")
            paths["fits"].write_bytes(b"previous new map")
            staged = directory / "staged"
            staged.write_bytes(b"replacement")
            with self.assertRaises(FileExistsError):
                producer.atomic_publish(staged, paths["fits"])
            self.assertEqual(paths["fits"].read_bytes(), b"previous new map")
            backup = producer.backup_existing(paths)
            producer.atomic_publish(staged, paths["fits"], overwrite=True)
            self.assertEqual(paths["fits"].read_bytes(), b"replacement")
            self.assertEqual(Path(backup["files"][0]["backup"]["path"]).read_bytes(), b"previous new map")
            self.assertEqual(paths["old_cmb"].read_bytes(), b"old map must remain")
            with producer.exclusive_task_lock(paths["lock"]):
                with self.assertRaises(FileExistsError):
                    with producer.exclusive_task_lock(paths["lock"]):
                        self.fail("Concurrent lock should fail")
                self.assertTrue(paths["lock"].exists())
            self.assertFalse(paths["lock"].exists())

    def test_partial_and_mismatched_existing_outputs_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp) / "00000"
            directory.mkdir()
            paths = producer.output_paths(temp, 0)
            self.assertIsNone(producer.verify_existing(paths, {}))
            paths["fits"].write_bytes(b"partial")
            with self.assertRaises(FileExistsError):
                producer.verify_existing(paths, {})
            paths["spectra"].write_bytes(b"spectrum")
            paths["json"].write_text("{}")
            with self.assertRaises(FileExistsError):
                producer.verify_existing(paths, {})
            paths["json"].write_text("invalid json")
            with self.assertRaises(FileExistsError):
                producer.verify_existing(paths, {})

    def test_temporary_fits_serialization_checksum_and_exact_resume(self):
        """IO-only nside8 fixture; this is not a production remapping test."""
        from astropy.io import fits
        with tempfile.TemporaryDirectory() as temp, patch.object(producer, "NSIDE", 8):
            directory = Path(temp) / "00000"
            directory.mkdir()
            paths = producer.output_paths(temp, 0)
            maps = np.random.default_rng(145).normal(size=(3, 12 * 8 ** 2))
            hp.write_map(paths["fits"], maps, dtype=np.float64, nest=False, coord="G",
                         fits_IDL=True, column_names=["T", "Q", "U"],
                         column_units=["uK_CMB"] * 3, extra_header=[("POLCCONV", "COSMO")])
            with fits.open(paths["fits"], mode="update") as hdus:
                for hdu in hdus:
                    hdu.add_checksum()
            serialized = producer._read_output_for_verification(paths["fits"], "G")
            self.assertEqual(serialized, [producer.full_statistics(field) for field in maps])
            np.savez_compressed(paths["spectra"], fixture=np.arange(5))
            contract = {"seed": 12, "coordinate_provenance": {"output_coordinate": "G"}}
            metadata = {"status": "COMPLETE", "contract": contract,
                        "contract_sha256": producer.canonical_sha256(contract),
                        "output": producer.file_record(paths["fits"]),
                        "input_spectra": producer.file_record(paths["spectra"]),
                        "finite_statistics": {"output_serialized_TQU": serialized}}
            paths["json"].write_text(json.dumps(metadata))
            self.assertTrue(producer.verify_existing(paths, contract)["reused"])
            with self.assertRaises(FileExistsError):
                producer.verify_existing(paths, {**contract, "seed": 13})
            with fits.open(paths["fits"], mode="update") as hdus:
                hdus[1].header["COORDSYS"] = "C"
            with self.assertRaises(FileExistsError):
                producer.verify_existing(paths, contract)

    def test_source_header_rejects_wrong_units_and_resolution(self):
        from astropy.io import fits
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "kappa.fits"
            table = fits.BinTableHDU.from_columns([fits.Column(name="T", format="D", array=np.zeros(3))])
            table.header["NSIDE"] = 2048
            table.header["ORDERING"] = "RING"
            table.header["LASTPIX"] = 12 * 2048 ** 2 - 1
            table.header["TUNIT1"] = "uK_CMB"
            fits.HDUList([fits.PrimaryHDU(), table]).writeto(path)
            with self.assertRaisesRegex(ValueError, "dimensionless"):
                producer.read_kappa_header(path)
            with fits.open(path, mode="update") as hdus:
                hdus[1].header["TUNIT1"] = "1"
            self.assertEqual(producer.read_kappa_header(path)["NSIDE"], 2048)
            with fits.open(path, mode="update") as hdus:
                hdus[1].header["NSIDE"] = 512
            with self.assertRaisesRegex(ValueError, "NSIDE2048"):
                producer.read_kappa_header(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
