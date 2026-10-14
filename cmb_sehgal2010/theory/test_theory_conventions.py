"""Cheap analytic checks of zero-mode export; no CAMB solve or sky transform."""
import importlib.util
import json
from pathlib import Path
import unittest

import numpy as np


def load_module(filename):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MODULES = [load_module('build_sehgal2010_theory.py'), load_module('relens_primary_with_pp.py')]


class ZeroModeConventionTests(unittest.TestCase):
    def test_nonzero_library_zero_mode_preserved_without_changing_positive_ell(self):
        ell = np.arange(9, dtype=np.int64)
        returned = np.arange(36, dtype=np.float64).reshape(9, 4) / 7
        returned[0] = [.0008265324, -.25, .125, 9.]
        before = returned.copy()
        factor = ell * (ell + 1) / (2 * np.pi)
        for module in MODULES:
            with self.subTest(module=module.__name__):
                raw, canonical_D, zero = module.canonicalize_correlation_spectra(np, returned, ell)
                np.testing.assert_array_equal(returned, before)
                np.testing.assert_array_equal(canonical_D[1:], returned[1:])
                np.testing.assert_array_equal(zero, returned[0])
                np.testing.assert_array_equal(raw[0], np.zeros(4))
                np.testing.assert_array_equal(canonical_D[0], np.zeros(4))
                np.testing.assert_allclose(raw * factor[:, None], canonical_D, rtol=1e-14, atol=0)
                typed_metadata = json.loads(json.dumps({'raw_correlation_zero_mode': zero.tolist(),
                                                       'canonical_raw_C0_and_D0': 0.0}, allow_nan=False))
                self.assertTrue(all(isinstance(value, float) for value in typed_metadata['raw_correlation_zero_mode']))

    def test_zero_lensing_primary_identity(self):
        ell = np.arange(9, dtype=np.int64)
        primary_raw = np.arange(36, dtype=np.float64).reshape(9, 4) / 100
        primary_raw[0] = 0
        primary_D = primary_raw * (ell * (ell + 1) / (2 * np.pi))[:, None]
        for module in MODULES:
            with self.subTest(module=module.__name__):
                _, canonical_D, zero = module.canonicalize_correlation_spectra(np, primary_D, ell)
                np.testing.assert_array_equal(canonical_D, primary_D)
                np.testing.assert_array_equal(zero, np.zeros(4))

    def test_nonfinite_or_wrong_shape_rejected(self):
        ell = np.arange(9, dtype=np.int64)
        invalid = np.zeros((9, 4))
        invalid[0, 0] = np.nan
        for module in MODULES:
            with self.subTest(module=module.__name__):
                with self.assertRaises(ValueError):
                    module.canonicalize_correlation_spectra(np, invalid, ell)
                with self.assertRaises(ValueError):
                    module.canonicalize_correlation_spectra(np, np.zeros((8, 4)), ell)


if __name__ == '__main__':
    unittest.main()
