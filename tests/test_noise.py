"""Small native-grid controls for the released scalar noise generator."""
from dataclasses import asdict, replace
import json
import tempfile
from pathlib import Path
import unittest
from unittest import mock

import healpy as hp
import numpy as np
from astropy.io import fits

from ali_sky import noise
from ali_sky.config import normalize
from ali_sky.io import read_fits


class NoiseTests(unittest.TestCase):
    def setUp(self):
        self.channel = noise.default_channels("LAT")[4]

    @staticmethod
    def legacy_draw(channel, nside, scale_factor, seeds):
        """Independent transcription of the archived obs_noise scalar law."""
        legacy_lmax = 3 * (nside * scale_factor) - 1
        ell = np.arange(1, legacy_lmax + 1)
        white = (np.pi / (60 * 180) * channel.nlev_uK_arcmin) ** 2
        red_p = white if channel.n_red_P is None else (np.pi / (60 * 180) * channel.n_red_P) ** 2
        if not channel.include_red:
            red_p, red_t, knee_t, knee_p, alpha_t, alpha_p = 0, 0, 1, 1, 1, 1
        else:
            red_t, knee_t, knee_p = channel.N_red_T, channel.l_knee_T, channel.l_knee_P
            alpha_t, alpha_p = channel.alpha_knee_T, channel.alpha_knee_P
        nt = red_t * (ell / knee_t) ** alpha_t + white
        npol = 2 * (red_p * (ell / knee_p) ** alpha_p + white)
        state = np.random.get_state()
        try:
            output = []
            for label, spectrum in zip(("T", "Q", "U"), (nt, npol, npol)):
                np.random.seed(seeds[label])
                output.append(hp.ud_grade(hp.synfast(spectrum, nside * scale_factor,
                                                   pol=False, new=True), nside))
            return np.array(output), nt, npol
        finally:
            np.random.set_state(state)

    def test_active_defaults(self):
        sat, lat = noise.default_channels("sat"), noise.default_channels("LAT")
        self.assertEqual([c.frequency_GHz for c in sat], [27, 40, 90, 150, 225, 280])
        self.assertEqual([c.nlev_uK_arcmin for c in sat], [18.27, 15.51, 3.60, 4.57, 11.46, 22.49])
        self.assertEqual([c.beam_fwhm_arcmin for c in sat], [97.34, 65.70, 16.20, 9.70, 10.60, 8.70])
        self.assertTrue(all(c.N_red_T == 0 for c in sat))
        self.assertEqual([c.l_knee_P for c in sat], [30, 30, 50, 50, 70, 100])
        self.assertEqual([c.alpha_knee_P for c in sat], [-2.4, -2.4, -2.5, -3, -3, -3])
        self.assertEqual([c.nlev_uK_arcmin for c in lat], [52, 27, 5.8, 6.3, 15, 37])
        self.assertEqual([c.beam_fwhm_arcmin for c in lat], [7.4, 5.1, 2.2, 1.4, 1, .9])
        self.assertEqual([c.N_red_T for c in lat], [v / (8760 * 3600) for v in [100, 39, 230, 1500, 17000, 31000]])
        self.assertTrue(all(c.l_knee_T == 1000 and c.alpha_knee_T == -3.5 for c in sat + lat))
        self.assertTrue(all(c.l_knee_P == 700 and c.alpha_knee_P == -1.4 for c in lat))

    def test_legacy_law_and_synfast_exact(self):
        for channel in (self.channel, noise.default_channels("SAT")[2],
                        replace(self.channel, n_red_P=None), replace(self.channel, include_red=False)):
            for scale_factor in (1, 2):
                with self.subTest(channel=channel, scale_factor=scale_factor):
                    maps, metadata = noise.generate_noise(channel, 4, 37, 20261011, "LAT", scale_factor)
                    expected, nt, npol = self.legacy_draw(channel, 4, scale_factor, metadata["seeds"])
                    actual_t, actual_p = noise.noise_spectra(channel, 4, scale_factor)
                    np.testing.assert_array_equal(actual_t, nt)
                    np.testing.assert_array_equal(actual_p, npol)
                    np.testing.assert_array_equal(maps, expected)
                    self.assertEqual(maps.shape, (3, hp.nside2npix(4)))
                    self.assertEqual(metadata["lmax"], 3 * 4 * scale_factor - 2)
                    self.assertEqual(metadata["spectrum_array_length"], 3 * 4 * scale_factor - 1)
                    self.assertEqual(metadata["statistics"][0]["finite_count"], hp.nside2npix(4))
                    json.dumps(metadata, allow_nan=False)

    def test_ell_plus_one_and_white_only(self):
        nt, npol = noise.noise_spectra(self.channel, 4)
        white = (np.pi / 10800 * 15) ** 2
        self.assertEqual(nt[0], self.channel.N_red_T * (1 / 1000) ** -3.5 + white)
        self.assertEqual(npol[0], 2 * (white * (1 / 700) ** -1.4 + white))
        white_t, white_p = noise.noise_spectra(replace(self.channel, include_red=False), 4)
        np.testing.assert_array_equal(white_t, np.full(11, white))
        np.testing.assert_array_equal(white_p, np.full(11, 2 * white))

    def test_seed_order_rank_and_id_independence(self):
        channels = noise.default_channels("LAT")[3:5]
        tasks = [(frequency, rid) for frequency in (0, 1) for rid in (37, 42)]
        baseline = {}
        for index, rid in tasks:
            c = channels[index]
            baseline[(c.frequency_GHz, rid)] = noise.generate_noise(c, 4, rid, 9876, "LAT")
        # A reversed channel list and two noncontiguous MPI-style partitions.
        reversed_channels = channels[::-1]
        partitioned = [(1, 42), (0, 37), (1, 37), (0, 42)]
        for index, rid in partitioned:
            c = reversed_channels[index]
            value, meta = noise.generate_noise(c, 4, rid, 9876, "lat")
            old_value, old_meta = baseline[(c.frequency_GHz, rid)]
            np.testing.assert_array_equal(value, old_value)
            self.assertEqual(meta["seeds"], old_meta["seeds"])
        self.assertFalse(np.array_equal(baseline[(150, 37)][0], baseline[(150, 42)][0]))
        self.assertNotEqual(baseline[(150, 37)][1]["seeds"], baseline[(225, 37)][1]["seeds"])
        old_seed = baseline[(225, 37)][1]["seeds"]
        paired = noise.generate_noise(replace(channels[1], nlev_uK_arcmin=25), 4, 37, 9876, "LAT")[1]
        self.assertEqual(old_seed, paired["seeds"])

    def test_versioned_seed_policy(self):
        _, metadata = noise.generate_noise(self.channel, 4, 37, 20261011, "LAT")
        self.assertEqual(metadata["seed_policy"], "ali-sky-noise-seedsequence-v1")
        self.assertEqual(metadata["seeds"], {"T": 2321517054, "Q": 2321518288, "U": 2321521375})
        for rid, master_seed in [(0, 0), (2 ** 64 - 1, 2 ** 64 - 1)]:
            _, meta = noise.generate_noise(self.channel, 4, rid, master_seed)
            self.assertEqual(meta["seeds"]["Q"], (meta["seeds"]["T"] + 1234) % 2 ** 32)
            self.assertEqual(meta["seeds"]["U"], (meta["seeds"]["T"] + 4321) % 2 ** 32)

    def test_rng_restored_success_and_failure(self):
        np.random.seed(445566)
        expected_state = np.random.get_state()
        noise.generate_noise(self.channel, 4, 0, 7654, "LAT")
        self.assertEqual(expected_state[0], np.random.get_state()[0])
        np.testing.assert_array_equal(expected_state[1], np.random.get_state()[1])
        self.assertEqual(expected_state[2:], np.random.get_state()[2:])
        with mock.patch.object(noise.hp, "synfast", side_effect=RuntimeError("test transform failure")):
            with self.assertRaisesRegex(RuntimeError, "test transform failure"):
                noise.generate_noise(self.channel, 4, 0, 7654, "LAT")
        np.testing.assert_array_equal(expected_state[1], np.random.get_state()[1])
        self.assertEqual(expected_state[2:], np.random.get_state()[2:])

    def test_invalid_parameters(self):
        for field, values in {"frequency_GHz": [0, -1, np.nan, np.inf, True],
                              "nlev_uK_arcmin": [-1, np.inf, False],
                              "beam_fwhm_arcmin": [-1], "N_red_T": [-1, np.nan],
                              "l_knee_T": [0], "l_knee_P": [0],
                              "alpha_knee_P": [np.inf], "include_red": [1]}.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(self.channel, **{field: value})
        with self.assertRaises(ValueError):
            noise.channel_from_dict({**asdict(self.channel), "typo": 5})
        with self.assertRaises(ValueError):
            noise.channel_from_dict({"frequency_GHz": 225})
        for nside, factor in [(0, 1), (3, 1), (4, 0), (4, 3), (True, 1), (4, 1.0)]:
            with self.subTest(nside=nside, factor=factor), self.assertRaises(ValueError):
                noise.noise_spectra(self.channel, nside, factor)
        for key, value in [("rid", -1), ("master_seed", -1), ("rid", True), ("master_seed", 2 ** 64)]:
            kwargs = {"channel": self.channel, "nside": 4, "rid": 0, "master_seed": 7, key: value}
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                noise.generate_noise(**kwargs)
        with self.assertRaises(ValueError):
            noise.noise_spectra(replace(self.channel, alpha_knee_T=-1e308), 4)
        with self.assertRaises(ValueError):
            noise.noise_spectra(replace(self.channel, nlev_uK_arcmin=1e308), 4)
        with self.assertRaises(ValueError):
            noise.generate_noise(replace(self.channel, nlev_uK_arcmin=1e308), 4, 0, 0)

    def test_zero_channel(self):
        zero = replace(self.channel, nlev_uK_arcmin=0, N_red_T=0, n_red_P=0)
        maps, meta = noise.generate_noise(zero, 4, 0, 0)
        np.testing.assert_array_equal(maps, np.zeros((3, 192)))
        self.assertTrue(all(s["min"] == s["max"] == 0 for s in meta["statistics"]))

    def test_fits_roundtrip_and_safe_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            config = normalize({"nside": 4, "seed": 9123, "frequencies_GHz": [225],
                                "start_id": 37, "n_realizations": 2, "output_dir": directory}, "noise")
            self.assertEqual(config["coordinate"], "G")
            written = noise.generate_noisefile(config, (0, 37))
            path = Path(directory) / "LAT" / "225" / "noise_00037.fits"
            self.assertEqual(written["output"]["path"], str(path))
            self.assertFalse(written["reused"])
            actual, coord, info = read_fits(path, 4)
            expected, meta = noise.generate_noise(self.channel, 4, 37, 9123, "LAT")
            np.testing.assert_array_equal(actual, expected)
            self.assertEqual(coord, "G")
            self.assertEqual(info["sha256"], written["output"]["sha256"])
            header = fits.getheader(path, 1)
            self.assertEqual(header["COORDSYS"], "G")
            self.assertEqual(header["ORDERING"], "RING")
            self.assertEqual(header["POLCCONV"], "COSMO")
            self.assertEqual([header[f"TTYPE{i}"] for i in range(1, 4)], ["T", "Q", "U"])
            self.assertEqual([header[f"TUNIT{i}"] for i in range(1, 4)], ["uK_CMB"] * 3)
            with fits.open(path, checksum=True) as hdul:
                self.assertTrue(all(hdu.verify_checksum() == 1 for hdu in hdul))
            with mock.patch.object(noise, "generate_noise", side_effect=AssertionError("reused task must not generate")):
                resumed = noise.generate_noisefile(config, (0, 37))
            self.assertTrue(resumed["reused"])
            self.assertEqual(resumed["seeds"], meta["seeds"])
            # These existing-file changes must be rejected, never silently reused.
            for key, value in [("seed", 9124), ("coordinate", "C"), ("nside", 8), ("scale_factor", 2)]:
                with self.subTest(key=key), self.assertRaises(FileExistsError):
                    noise.generate_noisefile({**config, key: value}, (0, 37))
            changed = dict(config)
            changed["channels"] = [{**config["channels"][0], "beam_fwhm_arcmin": 2}]
            with self.assertRaises(FileExistsError):
                noise.generate_noisefile(changed, (0, 37))
            original = path.with_suffix(".json").read_text()
            bad = json.loads(original)
            bad["seeds"]["Q"] += 1
            path.with_suffix(".json").write_text(json.dumps(bad))
            with self.assertRaises(FileExistsError):
                noise.generate_noisefile(config, (0, 37))
            path.with_suffix(".json").write_text(original)
            bad = json.loads(original)
            bad["output"] = []
            path.with_suffix(".json").write_text(json.dumps(bad))
            with self.assertRaises(FileExistsError):
                noise.generate_noisefile(config, (0, 37))
            path.with_suffix(".json").write_text(original)
            # Byte corruption must fail even if the numerical FITS map still reads.
            with path.open("r+b") as handle:
                handle.seek(-1, 2)
                last = handle.read(1)
                handle.seek(-1, 2)
                handle.write(bytes([last[0] ^ 1]))
            with self.assertRaises(FileExistsError):
                noise.generate_noisefile(config, (0, 37))
            replaced = noise.generate_noisefile({**config, "overwrite": True}, (0, 37))
            self.assertFalse(replaced["reused"])
            self.assertTrue(noise.generate_noisefile(config, (0, 37))["reused"])
            with self.assertRaises(ValueError):
                noise.generate_noisefile(config, (0, 99))


if __name__ == "__main__":
    unittest.main()
