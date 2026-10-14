"""Independent small native-grid controls for the portable forward model.

Every sky input below is a fresh NSIDE=8/16 fixture; production maps and
manuscript artifacts are never read, resampled or modified by this suite.
"""
from contextlib import contextmanager, redirect_stdout
import copy
import io as textio
import json
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

import healpy as hp
import numpy as np
from astropy.io import fits

from ali_sky import cli, config, coordinates, foregrounds, observation, sed
from ali_sky import io


def fixture_fits(path, maps, coordinate="G", unit="uK_CMB", pol="COSMO"):
    """Write fixture inputs directly with healpy, independently of ali_sky IO."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 1 if np.asarray(maps).ndim == 1 else 3
    hp.write_map(str(path), maps, dtype=np.float64, nest=False,
                 coord=coordinate, column_units=[unit] * count,
                 extra_header=[("POLCCONV", pol)], overwrite=True)
    return path


def independent_sed(component, frequency):
    """Transcribe archived physical constants, rather than call package SED."""
    nu = np.float64(frequency)
    if component == "ksz":
        return 1.
    if component in ("tsz", "radio"):
        x = 6.62607004e-34 * (nu * 1.e9) / (1.38064852e-23 * 2.726)
        if component == "tsz":
            return x * (np.exp(x) + 1.) / np.expm1(x) - 4.
        logx = np.log(x)
        return np.exp(-2.7714227059838192 * logx + .019242069853106277 * logx**2
                      + .36813523402646997 * logx**3)
    x = 6.626068e-34 * (nu * 1.e9) / (1.3806503e-23 * 2.726)
    derivative = (nu * 1.e9)**3 * np.exp(x) / np.expm1(x)**2 * x / 2.726
    dust_x = 6.626068e-34 * (nu * 1.e9) / (1.3806503e-23 * 8.35)
    return nu**(3. + 1.59) / np.expm1(dust_x) / derivative


class ObservationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.nside = 8
        self.npix = hp.nside2npix(self.nside)
        theta, phi = hp.pix2ang(self.nside, np.arange(self.npix))
        self.pattern = np.array([2. + np.cos(theta), .3 * np.sin(theta) * np.cos(phi),
                                 .2 * np.sin(theta) * np.sin(phi)])
        self.rid, self.noise_id = 7, 3
        self.sky_dir = self.root / "repository" / "sky"
        self.gal_dir = self.root / "repository" / "galactic"
        self.noise_dir = self.root / "provided_noise"
        self.mask_path = self.root / "mask.fits"
        self.mask = np.full(self.npix, .5)
        self.mask[:31] = 0.
        fixture_fits(self.mask_path, self.mask, unit="")
        self.cmb = self.pattern.copy()
        self.noise = self.pattern[:, ::-1].copy() * .07
        self.gal = {}
        self.xfg = {}
        folder = self.sky_dir / f"{self.rid:05d}"
        fixture_fits(folder / f"lensedcmb_new_{self.rid:05d}.fits", self.cmb)
        for index, (name, preset) in enumerate(config.PRESETS.items(), start=1):
            maps = self.pattern * (.01 * index)
            if name in ("ame", "freefree"):
                maps[1:] = 0.
            self.gal[name] = maps
            fixture_fits(self.gal_dir / f"fg_{name}_{preset}_150GHz_nside8_G.fits", maps)
        for index, name in enumerate(("tsz", "cib", "radio", "ksz"), start=1):
            maps = np.sin(theta) * (.03 * index) + .1 * index
            self.xfg[name] = maps
            filename = f"ksz_{self.rid:05d}.fits" if name == "ksz" else f"{observation.PREFIX[name]}_148ghz_{self.rid:05d}.fits"
            fixture_fits(folder / filename, maps)
        fixture_fits(self.noise_dir / f"noise_{self.noise_id:05d}.fits", self.noise)
        self.raw = {
            "instrument": "LAT", "nside": 8, "ids": [self.rid],
            "noise_ids": [self.noise_id], "coordinate": "G",
            "channels": [{"frequency_GHz": 150, "beam_fwhm_arcmin": 0}],
            "repository": {"sky_dir": str(self.sky_dir), "galactic_dir": str(self.gal_dir),
                           "noise_dir": str(self.root / "legacy_noise"), "mask": str(self.mask_path),
                           "sky_coordinate": "G", "noise_coordinate": "G", "mask_coordinate": "G"},
            "noise_dir": str(self.noise_dir), "noise_pattern": "noise_{id:05d}.fits",
            "mask": str(self.mask_path), "mask_coordinate": "G", "use_repository_mask": False,
            "cache_dir": str(self.root / "cache"), "output_dir": str(self.root / "output"),
            "lmax": 15, "harmonic_iter": 3,
        }

    def normalized(self, **changes):
        return config.normalize({**copy.deepcopy(self.raw), **changes}, "observe", self.root)

    def expected_sky(self, frequency=150):
        result = self.cmb + sum(self.gal.values())
        for name, maps in self.xfg.items():
            factor = independent_sed(name, frequency) / independent_sed(name, 148)
            result[0] += factor * maps
        return result

    def test_configuration_exact_frequencies_and_native_geometry(self):
        normal = self.normalized()
        self.assertEqual(normal["coordinate"], "G")
        self.assertIn("Dataset-owner declaration", normal["repository"]["frame_basis"])
        self.assertFalse(normal["custom_noise"])  # A beam-only override does not generate fresh noise.
        for bad in ({"nside": 12}, {"lmax": 24}, {"harmonic_iter": -1},
                    {"noise_ids": []}, {"ids": [7, 7]}, {"coordinate": "Z"}):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.normalized(**bad)
        for frequency in (39, 93, 145, 148):
            with self.subTest(frequency=frequency), self.assertRaisesRegex(ValueError, "explicit"):
                self.normalized(channels=None, frequencies_GHz=[frequency])
        custom = self.normalized(channels=[{"frequency_GHz": 148, "nlev_uK_arcmin": 12,
                                             "beam_fwhm_arcmin": 2}])
        self.assertEqual(custom["channels"][0]["frequency_GHz"], 148)
        self.assertTrue(custom["custom_noise"])

    def test_sed_independent_amplitude_ratios(self):
        for component in ("tsz", "cib", "radio", "ksz"):
            for frequency in (27, 40, 90, 148, 150, 225, 280):
                with self.subTest(component=component, frequency=frequency):
                    expected = independent_sed(component, frequency) / independent_sed(component, 148)
                    self.assertAlmostEqual(sed.scale_factor(component, frequency), expected,
                                           delta=abs(expected) * 2.e-13 + 1.e-15)
            self.assertEqual(sed.scale_factor(component, 148), 1.)
        self.assertLess(sed.scale_factor("tsz", 225), 0.)
        for frequency in (0, -1, np.inf, np.nan):
            with self.assertRaises(ValueError):
                sed.scale_factor("radio", frequency)

    def test_exact_archived_band_has_precedence_over_148(self):
        cfg = self.normalized()
        folder = self.sky_dir / f"{self.rid:05d}"
        for component in ("tsz", "cib", "radio"):
            native = self.xfg[component] + 9.
            fixture_fits(folder / f"{observation.PREFIX[component]}_150ghz_{self.rid:05d}.fits", native)
            maps, meta = observation.read_sky(cfg, component, self.rid, 150)
            np.testing.assert_array_equal(maps, native)
            self.assertEqual(meta["SED_factor"], 1.)
            self.assertEqual(meta["source_frequency_GHz"], 150)
        cmb, _ = observation.read_sky(cfg, "cmb", self.rid)
        np.testing.assert_array_equal(cmb, self.cmb)

    def test_nonwarehouse_148_uses_true_source_and_custom_galactic_band(self):
        cfg = self.normalized(channels=[{"frequency_GHz": 148, "nlev_uK_arcmin": 12,
                                        "beam_fwhm_arcmin": 0}], add_noise=False)
        for name, preset in config.PRESETS.items():
            fixture_fits(self.gal_dir / f"fg_{name}_{preset}_148GHz_nside8_G.fits", self.gal[name])
        for name in ("tsz", "cib", "radio"):
            maps, meta = observation.read_sky(cfg, name, self.rid, 148)
            np.testing.assert_array_equal(maps, self.xfg[name])
            self.assertEqual(meta["SED_factor"], 1.)
        result = observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)
        maps = hp.read_map(result["output"]["path"], field=(0, 1, 2), dtype=np.float64)
        np.testing.assert_allclose(maps, self.expected_sky(148) * self.mask, rtol=2.e-14, atol=2.e-14)

    def test_sehgal_release_cmb_ignores_discrepant_legacy_map(self):
        cfg = self.normalized()
        folder = self.sky_dir / f"{self.rid:05d}"
        legacy = fixture_fits(folder / f"lensedcmb_{self.rid:05d}.fits", self.cmb + 1234.)
        current = folder / f"lensedcmb_new_{self.rid:05d}.fits"
        observation.preflight(cfg)
        maps, metadata = observation.read_sky(cfg, "cmb", self.rid)
        np.testing.assert_array_equal(maps, self.cmb)
        self.assertEqual(metadata["path"], str(current.resolve()))
        self.assertEqual(metadata["sha256"], io.sha256(current))
        self.assertNotEqual(metadata["sha256"], io.sha256(legacy))
        result = observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)
        actual = hp.read_map(result["output"]["path"], field=(0, 1, 2), dtype=np.float64)
        np.testing.assert_allclose(actual, (self.expected_sky() + self.noise) * self.mask,
                                   rtol=2.e-14, atol=2.e-14)
        self.assertEqual(result["inputs"]["cmb"]["path"], str(current.resolve()))

    def test_legacy_cmb_does_not_replace_missing_sehgal_release_map(self):
        cfg = self.normalized()
        folder = self.sky_dir / f"{self.rid:05d}"
        fixture_fits(folder / f"lensedcmb_{self.rid:05d}.fits", self.cmb + 1234.)
        (folder / f"lensedcmb_new_{self.rid:05d}.fits").unlink()
        with self.assertRaisesRegex(FileNotFoundError, "lensedcmb_new_00007"):
            observation.preflight(cfg)
        with self.assertRaisesRegex(FileNotFoundError, "lensedcmb_new_00007"):
            observation.read_sky(cfg, "cmb", self.rid)
        self.assertFalse(Path(cfg["output_dir"]).exists())

    def test_forward_sum_matched_noise_and_mask_applied_once(self):
        cfg = self.normalized()
        observation.preflight(cfg)
        result = observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)
        maps = hp.read_map(result["output"]["path"], field=(0, 1, 2), dtype=np.float64)
        np.testing.assert_allclose(maps, (self.expected_sky() + self.noise) * self.mask,
                                   rtol=2.e-14, atol=2.e-14)
        self.assertEqual(result["noise_id"], 3)
        self.assertEqual(result["sky_id"], 7)
        self.assertTrue(result["no_downgrade"])
        self.assertFalse(result["inputs"]["noise"]["beam_applied"])
        self.assertEqual(result["unit"], "uK_CMB")
        self.assertEqual(result["coordinate"], "G")
        header = fits.getheader(result["output"]["path"], 1)
        self.assertEqual(header["FREQGHZ"], 150.)
        self.assertEqual(header["ORDERING"], "RING")
        self.assertEqual(header["POLCCONV"], "COSMO")
        self.assertEqual(header["CREATOR"], "ali_sky 1.1.0")
        with self.assertRaises(FileExistsError):
            observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)

    def test_missing_noise_id_not_wrapped_or_substituted(self):
        cfg = self.normalized(noise_ids=[53])
        with self.assertRaisesRegex(FileNotFoundError, "noise_00053"):
            observation.preflight(cfg)
        self.assertFalse(Path(cfg["output_dir"]).exists())

    def test_explicit_output_c_aligns_every_input_and_mask(self):
        cfg = self.normalized(coordinate="C")
        result = observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)
        rotation = hp.Rotator(coord=["G", "C"])
        sky_and_noise = np.asarray(rotation.rotate_map_pixel(self.expected_sky() + self.noise))
        source_vectors = hp.Rotator(coord=["C", "G"])(np.asarray(hp.pix2vec(8, np.arange(self.npix))))
        expected = sky_and_noise * self.mask[hp.vec2pix(8, *source_vectors)]
        actual = hp.read_map(result["output"]["path"], field=(0, 1, 2), dtype=np.float64)
        np.testing.assert_allclose(actual, expected, rtol=2.e-14, atol=2.e-14)
        self.assertEqual(result["coordinate"], "C")
        for metadata in list(result["inputs"]["galactic"].values()) + [result["inputs"]["cmb"]] + list(result["inputs"]["extragalactic"].values()):
            self.assertEqual(metadata["source_coordinate"], "G")
            self.assertEqual(metadata["output_coordinate"], "C")
            self.assertEqual(metadata["coordinate_operation"], "spin-aware pixel interpolation")
        self.assertEqual(result["inputs"]["noise"]["source_coordinate"], "G")
        self.assertEqual(result["mask"]["source_coordinate"], "G")
        self.assertEqual(result["mask"]["output_coordinate"], "C")

    def test_native_input_resolution_and_invalid_pixels_rejected(self):
        cfg = self.normalized()
        path = self.sky_dir / f"{self.rid:05d}" / f"lensedcmb_new_{self.rid:05d}.fits"
        fixture_fits(path, np.zeros((3, hp.nside2npix(16))))
        with self.assertRaisesRegex(ValueError, "no automatic resampling"):
            observation.read_sky(cfg, "cmb", self.rid)
        for bad in (np.nan, hp.UNSEEN):
            value = self.cmb.copy()
            value[1, 5] = bad
            fixture_fits(path, value)
            with self.assertRaisesRegex(ValueError, "non-finite|UNSEEN"):
                observation.read_sky(cfg, "cmb", self.rid)

    def test_fits_frame_units_and_owner_declaration_provenance(self):
        cfg = self.normalized()
        path = self.sky_dir / f"{self.rid:05d}" / f"lensedcmb_new_{self.rid:05d}.fits"
        with fits.open(path, mode="update") as hdus:
            del hdus[1].header["COORDSYS"]
            for i in (1, 2, 3):
                del hdus[1].header[f"TUNIT{i}"]
        maps, meta = observation.read_sky(cfg, "cmb", self.rid)
        np.testing.assert_array_equal(maps, self.cmb)
        self.assertFalse(meta["source_frame_verified"])
        self.assertEqual(meta["coordinate_operation"], "none")
        self.assertEqual(meta["source_coordinate"], "G")
        self.assertIn("owner declaration", meta["frame_basis"])
        self.assertIn("ICRS", meta["historical_source_frame"])
        with self.assertRaisesRegex(ValueError, "explicit source frame"):
            io.read_fits(path, 8, assume_unit="uK_CMB")
        fixture_fits(path, self.cmb, "C")
        maps, meta = observation.read_sky(cfg, "cmb", self.rid)
        np.testing.assert_allclose(maps, hp.Rotator(coord=["C", "G"]).rotate_map_pixel(self.cmb),
                                   rtol=0, atol=1.e-14)
        self.assertTrue(meta["source_frame_verified"])
        self.assertEqual(meta["source_coordinate"], "C")
        for ambiguous_unit in ("uK_RJ", "uK", " uK "):
            fixture_fits(path, self.cmb, unit=ambiguous_unit)
            with self.subTest(unit=ambiguous_unit), self.assertRaisesRegex(ValueError, "thermodynamic"):
                observation.read_sky(cfg, "cmb", self.rid)
        fixture_fits(path, self.cmb, pol="IAU")
        with self.assertRaisesRegex(ValueError, "polarization convention"):
            observation.read_sky(cfg, "cmb", self.rid)

    def test_coordinate_transfer_matches_official_spin_basis(self):
        self.assertIs(coordinates.rotate_map(self.pattern, "G", "G"), self.pattern)
        for source, target in (("G", "C"), ("C", "G"), ("G", "E")):
            rotator = hp.Rotator(coord=[source, target])
            expected = np.asarray(rotator.rotate_map_pixel(self.pattern))
            actual = coordinates.rotate_map(self.pattern, source, target, chunk=101)
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1.e-14)
            scalar = coordinates.rotate_map(self.pattern[0], source, target, chunk=83)
            np.testing.assert_allclose(scalar, expected[0], rtol=0, atol=1.e-14)
            theta, phi = rotator.I(*hp.pix2ang(8, np.arange(self.npix)))
            wrong_q = hp.get_interp_val(self.pattern[1], theta, phi)
            self.assertGreater(float(np.max(np.abs(actual[1] - wrong_q))), .05)

    def test_mask_frame_explicit_and_nearest_values_preserved(self):
        with fits.open(self.mask_path, mode="update") as hdus:
            del hdus[1].header["COORDSYS"]
        with self.assertRaisesRegex(ValueError, "explicit source frame"):
            observation.read_mask(self.normalized(mask_coordinate=None))
        cfg = self.normalized(mask_coordinate="C")
        actual, metadata = observation.read_mask(cfg)
        target_vectors = np.asarray(hp.pix2vec(8, np.arange(self.npix)))
        source_vectors = hp.Rotator(coord=["G", "C"])(target_vectors)
        indices = hp.vec2pix(8, *source_vectors)
        np.testing.assert_array_equal(actual, self.mask[indices])
        self.assertEqual(metadata["coordinate_operation"], "nearest-pixel mask rotation")
        self.assertTrue(set(np.unique(actual)) <= {0., .5})
        invalid = self.mask.copy()
        invalid[2] = 1.01
        fixture_fits(self.mask_path, invalid, unit="")
        with self.assertRaisesRegex(ValueError, "within"):
            observation.preflight(self.normalized())

    def test_legacy_noise_full_ring_indices_and_coordinate_required(self):
        path = self.root / "legacy.npy"
        legacy = np.zeros((3, 2, self.npix))
        legacy[:, 0] = np.arange(self.npix)
        legacy[:, 1] = self.noise
        np.save(path, legacy)
        actual, frame, _ = io.read_noise(path, 8, assume_coord="G")
        np.testing.assert_array_equal(actual, self.noise)
        self.assertEqual(frame, "G")
        with self.assertRaisesRegex(ValueError, "coordinate"):
            io.read_noise(path, 8)
        for bad_index in (3., 3.5):
            broken = legacy.copy()
            broken[1, 0, 4] = bad_index
            np.save(path, broken)
            with self.assertRaisesRegex(ValueError, "complete monotone RING"):
                io.read_noise(path, 8, assume_coord="G")
        np.save(path, self.noise)
        actual, _, _ = io.read_noise(path, 8, assume_coord="G")
        np.testing.assert_array_equal(actual, self.noise)

    def test_joint_spin_beam_and_custom_teb_window(self):
        nside, lmax = 16, 12
        alms = np.zeros((3, hp.Alm.getsize(lmax)), complex)
        alms[0, hp.Alm.getidx(lmax, 3, 1)] = 1. + .3j
        alms[1, hp.Alm.getidx(lmax, 4, 2)] = .7 - .2j
        alms[2, hp.Alm.getidx(lmax, 5, 3)] = .05 + .03j
        maps = hp.alm2map(alms, nside, lmax=lmax, pol=True)
        fwhm = 1800.
        measured = hp.map2alm(maps, lmax=lmax, iter=3, pol=True)
        ell = np.arange(lmax + 1)
        sigma = fwhm * np.pi / 10800. / np.sqrt(8. * np.log(2.))
        window = np.exp(-.5 * ell * (ell + 1) * sigma**2)
        filtered = np.array([hp.almxfl(field, window) for field in measured])
        expected = hp.alm2map(filtered, nside, lmax=lmax, pol=True)
        result, metadata = observation.beam_sky(maps, {"beam_fwhm_arcmin": fwhm}, lmax, 3)
        np.testing.assert_allclose(result, expected, rtol=2.e-14, atol=2.e-14)
        self.assertEqual(metadata["applied_to"], "sky only")
        self.assertFalse(metadata["pixel_window"])
        custom = np.stack([window, .8 * window, .3 * window])
        path = self.root / "beam.npy"
        np.save(path, custom)
        expected = hp.alm2map([hp.almxfl(field, transfer) for field, transfer in zip(measured, custom)],
                             nside, lmax=lmax, pol=True)
        result, metadata = observation.beam_sky(maps, {"beam_transfer_path": str(path)}, lmax, 3)
        np.testing.assert_allclose(result, expected, rtol=2.e-14, atol=2.e-14)
        self.assertEqual(metadata["columns"], "T,E,B")
        for value in (np.ones((2, 13)), np.ones(3), np.full(13, np.nan)):
            np.save(path, value)
            with self.assertRaises(ValueError):
                observation.beam_sky(maps, {"beam_transfer_path": str(path)}, lmax, 3)

    def test_noise_not_beamed_and_mask_not_squared(self):
        cfg = self.normalized(channels=[{"frequency_GHz": 150, "beam_fwhm_arcmin": 1800}])
        with mock.patch.object(observation, "read_sky", side_effect=lambda c, name, rid, frequency=None:
                               (np.zeros_like(self.cmb) if name == "cmb" else np.zeros(self.npix), {})), \
             mock.patch.object(observation, "get_galactic_component", side_effect=lambda *args:
                               (np.zeros_like(self.cmb), "G", {})):
            result = observation.observe_one(cfg, cfg["channels"][0], self.rid, self.noise_id)
        actual = hp.read_map(result["output"]["path"], field=(0, 1, 2), dtype=np.float64)
        np.testing.assert_array_equal(actual, self.noise * self.mask)

    @contextmanager
    def fake_pysm(self, unpolarized=False):
        """A tiny template-backed emission model; no real template download."""
        fixture_dir = self.root / "pysm_templates"
        fixture_dir.mkdir(exist_ok=True)
        template = fixture_dir / "fixture_template.bin"
        template.write_bytes(b"known-template-input-v1")
        calls = []
        output = self.pattern.copy()
        module = types.ModuleType("pysm3")
        module.__path__ = []
        units = types.ModuleType("pysm3.units")
        units.GHz, units.uK_CMB = 1., "thermodynamic-fixture-unit"
        units.cmb_equivalencies = lambda frequency: ("cmb", frequency)
        utils = types.ModuleType("pysm3.utils")
        utils.__path__ = []
        data = types.ModuleType("pysm3.utils.data")

        class RemoteData:
            def get(self, logical_name):
                return str(Path(os.environ["PYSM_LOCAL_DATA"]) / logical_name)

        class Quantity:
            def to(self, unit, equivalencies):
                if unit != units.uK_CMB or equivalencies[0] != "cmb":
                    raise AssertionError("Thermodynamic conversion not requested")
                return types.SimpleNamespace(value=output[:1] if unpolarized else output)

        class Sky:
            def __init__(self, nside, preset_strings):
                calls.append((nside, preset_strings))
                RemoteData().get("fixture_template.bin")

            def get_emission(self, frequency):
                calls.append(frequency)
                return Quantity()

        data.RemoteData, module.Sky, module.units = RemoteData, Sky, units
        modules = {"pysm3": module, "pysm3.units": units, "pysm3.utils": utils,
                   "pysm3.utils.data": data}
        with mock.patch.dict(sys.modules, modules), \
             mock.patch.dict(os.environ, {"PYSM_LOCAL_DATA": str(fixture_dir)}), \
             mock.patch.object(foregrounds.importlib.metadata, "version", return_value="fixture-1"):
            yield template, calls, output

    def test_galactic_repository_exact_frequency_and_preset(self):
        cfg = self.normalized()
        with mock.patch.object(foregrounds, "_cache_key", side_effect=AssertionError("PySM unnecessary")):
            for component in config.PRESETS:
                actual, frame, metadata = foregrounds.get_galactic_component(cfg, component, 150)
                np.testing.assert_array_equal(actual, self.gal[component])
                self.assertEqual(frame, "G")
                self.assertEqual(metadata["preset"], config.PRESETS[component])
                self.assertEqual(metadata["method"], "repository component")

    def test_fresh_pysm_cache_hit_metadata_and_tamper_rejection(self):
        cfg = self.normalized()
        with self.fake_pysm() as (template, calls, expected):
            actual, frame, first = foregrounds.get_galactic_component(cfg, "dust", 148)
            np.testing.assert_array_equal(actual, expected)
            self.assertEqual(frame, "G")
            self.assertEqual(first["method"], "PySM generated")
            self.assertEqual(first["templates"][0]["sha256"], io.sha256(template))
            count = len(calls)
            again, frame, second = foregrounds.get_galactic_component(cfg, "dust", 148)
            np.testing.assert_array_equal(again, actual)
            self.assertEqual(second["method"], "PySM cache hit")
            self.assertEqual(len(calls), count)
            cache = Path(first["input"]["path"])
            sidecar = cache.with_suffix(".json")
            stored = json.loads(sidecar.read_text())
            self.assertEqual(stored["beam"], "none")
            self.assertEqual(stored["bandpass"], "monochromatic")
            self.assertEqual(stored["cache_key"]["frequency_GHz"], 148.)
            self.assertEqual(stored["cache_key"]["pysm3_version"], "fixture-1")
            self.assertEqual(fits.getheader(cache, 1)["COORDSYS"], "G")
            stored["output"]["sha256"] = "0" * 64
            sidecar.write_text(json.dumps(stored))
            with self.assertRaisesRegex(ValueError, "cache identity mismatch"):
                foregrounds.get_galactic_component(cfg, "dust", 148)

    def test_cache_template_identity_and_context_are_checked(self):
        cfg = self.normalized()
        with self.fake_pysm() as (template, calls, expected):
            _, _, first = foregrounds.get_galactic_component(cfg, "dust", 148)
            template.write_bytes(b"different-template-content")
            with self.assertRaisesRegex(ValueError, "template|Template"):
                foregrounds.get_galactic_component(cfg, "dust", 148)
            template.write_bytes(b"known-template-input-v1")
            other = self.root / "other_templates"
            other.mkdir()
            (other / template.name).write_bytes(b"other-template-tree")
            cfg["template_dir"] = str(other)
            _, _, changed = foregrounds.get_galactic_component(cfg, "dust", 148)
            self.assertNotEqual(first["input"]["path"], changed["input"]["path"])
            self.assertEqual(changed["method"], "PySM generated")

    def test_unpolarized_pysm_emission_has_explicit_zero_q_u(self):
        cfg = self.normalized()
        with self.fake_pysm(unpolarized=True) as (template, calls, expected):
            maps, frame, metadata = foregrounds.get_galactic_component(cfg, "ame", 148)
            self.assertEqual(maps.shape, (3, self.npix))
            np.testing.assert_array_equal(maps[0], expected[0])
            np.testing.assert_array_equal(maps[1:], np.zeros((2, self.npix)))
            self.assertEqual(metadata["preset"], "a1")
            self.assertEqual(frame, "G")
            cached, _, _ = foregrounds.get_galactic_component(cfg, "ame", 148)
            np.testing.assert_array_equal(cached, maps)

    def test_cli_automatic_noise_notice_ids_and_parallel_configuration(self):
        cfg = self.normalized(noise_dir=None, ids=[7, 8], noise_ids=[3, 3],
                              channels=[{"frequency_GHz": 150, "nlev_uK_arcmin": 9,
                                         "beam_fwhm_arcmin": 0}])
        printed = textio.StringIO()
        dummy = {"output": {"path": "fixture-product"}}
        with mock.patch.object(cli, "under_mpi", return_value=False), \
             mock.patch.object(observation, "preflight") as preflight, \
             mock.patch.object(cli, "launch_noise", return_value=0) as launcher, \
             mock.patch.object(observation, "observe_one", return_value=dummy) as worker, \
             redirect_stdout(printed):
            self.assertEqual(cli.run_observe(cfg), 0)
        self.assertIn("NOTICE", printed.getvalue())
        self.assertIn("5 MPI ranks", printed.getvalue())
        self.assertIn("beam is not applied to noise", printed.getvalue())
        preflight.assert_called_once_with(cfg)
        generated = launcher.call_args.args[0]
        self.assertEqual(generated["ids"], [3])
        self.assertEqual(generated["coordinate"], "G")
        self.assertEqual(generated["mpi_ranks"], 5)
        self.assertEqual([c.args[2:] for c in worker.call_args_list], [(7, 3), (8, 3)])
        self.assertTrue(all("generated_noise" in c.args[0]["noise_dir"] for c in worker.call_args_list))
        campaign = json.loads((Path(cfg["output_dir"]) / "observation_campaign.json").read_text())
        self.assertEqual(campaign["status"], "COMPLETE")

    def test_cli_existing_or_disabled_noise_not_regenerated_failure_stops(self):
        dummy = {"output": {"path": "fixture-product"}}
        for cfg in (self.normalized(), self.normalized(noise_dir=None, add_noise=False)):
            with mock.patch.object(cli, "under_mpi", return_value=False), \
                 mock.patch.object(observation, "preflight"), \
                 mock.patch.object(cli, "launch_noise", side_effect=AssertionError("unexpected generation")), \
                 mock.patch.object(observation, "observe_one", return_value=dummy):
                self.assertEqual(cli.run_observe(cfg), 0)
        cfg = self.normalized(noise_dir=None, channels=[{"frequency_GHz": 150, "nlev_uK_arcmin": 9}])
        with mock.patch.object(cli, "under_mpi", return_value=False), \
             mock.patch.object(observation, "preflight"), \
             mock.patch.object(cli, "launch_noise", return_value=7), \
             mock.patch.object(observation, "observe_one") as worker:
            self.assertEqual(cli.run_observe(cfg), 7)
            worker.assert_not_called()
        with mock.patch.object(cli, "under_mpi", return_value=True):
            with self.assertRaisesRegex(ValueError, "serial coordinator"):
                cli.run_observe(cfg)

    def test_cli_mpi_subprocess_uses_matching_runtime(self):
        cfg = config.normalize({"nside": 8, "ids": [3], "frequencies_GHz": [150],
                                "output_dir": str(self.root / "mpi_noise")}, "noise", self.root)
        launcher = self.root / "mpiexec"
        launcher.write_text("fixture launcher; never executed")
        with mock.patch.dict(os.environ, {"ALI_SKY_MPIEXEC": str(launcher)}), \
             mock.patch.object(cli, "under_mpi", return_value=False), \
             mock.patch.object(cli.subprocess, "run", return_value=types.SimpleNamespace(returncode=0)) as run:
            self.assertEqual(cli.launch_noise(cfg), 0)
        command = run.call_args.args[0]
        self.assertEqual(command[:6], [str(launcher), "-n", "5", sys.executable, "-m", "ali_sky.cli"])
        self.assertEqual(command[6:8], ["_noise-worker", "--config"])
        self.assertEqual(json.loads(Path(command[-1]).read_text()), cfg)


if __name__ == "__main__":
    unittest.main()
