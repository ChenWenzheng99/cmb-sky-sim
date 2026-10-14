"""Small analytic-array tests; these are not low-resolution sky validations."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.stats import t

from spectral_statistics import bin_native, evaluate, numerical_iteration_control
from validate_inputs import canonical_sha256, check_theory, positive_semidefinite, sha256


POLICY = json.loads(Path(__file__).with_name('validation_policy.json').read_text())


def analytic_reference(nell=1025):
    ell = np.arange(nell)
    shape = (ell + 1.)**-2
    return np.array([3 * shape, .6 * shape, .03 * shape,
                     .4 * np.sqrt(3 * .6) * shape * np.sin(ell / 35.)])


def balanced_pilot():
    reference = analytic_reference()
    fluctuation = np.linspace(-.005, .005, 10)
    raw = reference[None] * (1 + fluctuation[:, None, None])
    raw[:, 3] += fluctuation[:, None] * np.sqrt(reference[0] * reference[1])[None]
    return reference, raw


class StatisticsTests(unittest.TestCase):
    def test_weighted_bins_and_D_not_C(self):
        ell = np.arange(1025)
        actual = bin_native(ell**2)
        expected = sum((2*l + 1) * l*l for l in range(25, 50)) / sum(2*l + 1 for l in range(25, 50))
        self.assertEqual(actual.shape, (39,))
        self.assertEqual(actual[0], expected)
        self.assertNotEqual(actual[0], np.mean(ell[25:50]**2))
        cl = analytic_reference()
        factor = ell * (ell + 1.) / (2 * np.pi)
        self.assertFalse(np.allclose(bin_native(cl * factor), bin_native(cl) * bin_native(factor)))

    def test_signed_TE_and_predeclared_t9_limits(self):
        reference, raw = balanced_pilot()
        report, arrays = evaluate(raw, reference, POLICY)
        self.assertEqual(report['primary_spectral_gate'], 'PASS')
        self.assertFalse(report['metrics']['TE']['TE_ratio_used'])
        self.assertTrue(np.isfinite(arrays['TE_residual_over_sqrt_TT_EE']).all())
        self.assertEqual(report['amplitude_t_limit_df9'], t.ppf(1 - .005 / 6, 9))
        self.assertEqual(report['shape_t_limit_df9'], t.ppf(1 - .005 / 234, 9))
        self.assertFalse(report['BB_is_formal_pure_CAMB_acceptance_gate'])

    def test_primary_bias_and_duplicate_rejection(self):
        reference, raw = balanced_pilot()
        shifted = raw.copy()
        shifted[:, 0] *= 1.1
        report, _ = evaluate(shifted, reference, POLICY)
        self.assertEqual(report['primary_spectral_gate'], 'REVIEW_REQUIRED')
        with self.assertRaises(ValueError):
            evaluate(np.repeat(reference[None], 10, axis=0), reference, POLICY)

    def test_sample_SD_not_SEM_and_iteration(self):
        reference, raw = balanced_pilot()
        _, arrays = evaluate(raw, reference, POLICY)
        factor = np.arange(1025) * (np.arange(1025) + 1.) / (2 * np.pi)
        expected = bin_native(raw * factor).std(0, ddof=1)
        np.testing.assert_array_equal(arrays['sample_std_D_b'], expected)
        report = numerical_iteration_control(raw[0], raw[0] * (1 + 1e-9), arrays, POLICY)
        self.assertEqual(report['status'], 'PASS')
        report = numerical_iteration_control(raw[0], raw[0] * 1.01, arrays, POLICY)
        self.assertEqual(report['status'], 'REVIEW_REQUIRED')


class TheoryFixtureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        theory = self.root / 'theory'
        self.out = theory / 'outputs'
        sources = theory / 'source_evidence'
        self.out.mkdir(parents=True)
        sources.mkdir()
        config = json.loads(Path(__file__).resolve().parents[1].joinpath('theory/sehgal2010_config.json').read_text())
        config_path = theory / 'sehgal2010_config.json'
        config_path.write_text(json.dumps(config))
        source = sources / 'analytic_fixture.txt'
        source.write_text('Synthetic provenance fixture, not a cosmological calculation.\n')
        records = [{'filename': source.name, 'sha256': sha256(source)}]
        source_manifest = sources / 'source_manifest.json'
        source_manifest.write_text(json.dumps(records))
        ell = np.arange(6501, dtype=np.int64)
        factor = ell * (ell + 1.) / (2 * np.pi)
        raw = analytic_reference(6501).T
        unlensed = raw.copy()
        unlensed[:, 2] = 0
        producer = unlensed.copy()
        producer[6000:] = 0
        phi = 1e-9 / (ell + 1.)**4
        pp = phi * (ell * (ell + 1.))**2 / (2 * np.pi)
        cutpp = pp.copy()
        cutpp[6144:] = 0
        self.arrays = {'ell': ell, 'unlensed_raw_cl': unlensed, 'lensed_raw_cl': raw,
                       'unlensed_D_ell': unlensed * factor[:, None], 'lensed_D_ell': raw * factor[:, None],
                       'producer_unlensed_raw_cl': producer, 'phi_raw_cl': phi, 'phi_PP': pp,
                       'producer_cutoff_ell': np.arange(1025, dtype=np.int64), 'producer_cutoff_phi_PP': cutpp,
                       'producer_cutoff_lensed_raw_cl': raw[:1025], 'producer_cutoff_lensed_D_ell': raw[:1025] * factor[:1025, None],
                       'producer_cutoff_sampling2_raw_cl': raw[:1025], 'producer_cutoff_sampling2_D_ell': raw[:1025] * factor[:1025, None],
                       'producer_cutoff_zero_phi_D_ell': producer[:1025] * factor[:1025, None]}
        self.manifest = {'schema': 'sehgal2010_camb_theory_v1', 'status': 'COMPLETE_THEORY_WITH_DECLARED_ASSUMPTIONS',
                         'configuration': config, 'configuration_canonical_sha256': canonical_sha256(config),
                         'configuration_file': {'sha256': sha256(config_path)}, 'versions': {'camb': '1.6.1'},
                         'producer_cutoffs': config['producer_cutoffs'],
                         'array_contract': {key: {'shape': list(value.shape), 'dtype': str(value.dtype),
                                                 'units': 'dimensionless' if value.ndim == 1 else 'uK_CMB^2',
                                                 **({'columns': ['TT', 'EE', 'BB', 'TE']} if value.ndim == 2 else {})}
                                            for key, value in self.arrays.items()},
                         'normalization': {'trial_As': 2.1e-9, 'trial_sigma8': .8, 'As_calibrated': 2.1e-9,
                                           'achieved_sigma8_linear_z0': .8, 'relative_error': 0., 'full_solve_count': 2,
                                           'same_numerical_settings_in_both_solves': True, 'nonlinear_lensing_recalculated_at_final_As': True},
                         'source_manifest_sha256': sha256(source_manifest), 'source_records': records}
        self.write()

    def tearDown(self):
        self.temp.cleanup()

    def write(self):
        archive = self.out / 'sehgal2010_theory.npz'
        np.savez_compressed(archive, **self.arrays)
        self.manifest['artifact_manifest'] = {archive.name: {'sha256': sha256(archive), 'bytes': archive.stat().st_size}}
        (self.out / 'theory_manifest.json').write_text(json.dumps(self.manifest))

    def test_consistent_synthetic_contract(self):
        self.assertEqual(check_theory(self.root)['status'], 'PASS')

    def test_bad_column_contract_rejected(self):
        self.manifest['array_contract']['unlensed_raw_cl']['columns'] = ['TT', 'TE', 'EE', 'BB']
        self.write()
        with self.assertRaises(ValueError):
            check_theory(self.root)

    def test_unlensed_BB_corruption_rejected(self):
        self.arrays['unlensed_raw_cl'][100, 2] = 1e-5
        self.arrays['unlensed_D_ell'] = self.arrays['unlensed_raw_cl'] * (self.arrays['ell'] * (self.arrays['ell'] + 1.) / (2 * np.pi))[:, None]
        self.write()
        with self.assertRaises(ValueError):
            check_theory(self.root)

    def test_PP_normalization_corruption_rejected(self):
        self.arrays['phi_PP'] *= 2 * np.pi
        self.write()
        with self.assertRaises(AssertionError):
            check_theory(self.root)

    def test_raw_TE_covariance_rejected(self):
        invalid = analytic_reference().T
        invalid[:, 3] = 100
        with self.assertRaises(ValueError):
            positive_semidefinite(invalid, 'corrupted raw Cl')


if __name__ == '__main__':
    unittest.main()
