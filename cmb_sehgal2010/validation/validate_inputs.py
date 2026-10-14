#!/usr/bin/env python3
"""Independent, CAMB-free checks of the new theory and producer inputs.

This module reads small spectra and provenance files only. It never reads a
native map, computes a spherical transform, or calls CAMB/lenspyx.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


ORDER = ['TT', 'EE', 'BB', 'TE']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_record(path, record):
    path = Path(path)
    require(path.is_file(), f'Missing input: {path}')
    actual = {'bytes': path.stat().st_size, 'sha256': sha256(path)}
    require(actual['sha256'] == record['sha256'], f'Hash mismatch: {path}')
    if 'bytes' in record:
        require(actual['bytes'] == record['bytes'], f'Size mismatch: {path}')
    return {'path': str(path.resolve()), **actual}


def positive_semidefinite(values, name):
    require(np.isfinite(values).all(), f'Nonfinite {name}')
    require(np.all(values[:, :3] >= 0), f'Negative auto-spectrum: {name}')
    product = values[:, 0] * values[:, 1]
    require(np.all(values[:, 3]**2 <= product * (1 + 1e-10) + 1e-300),
            f'Invalid T/E covariance: {name}')


def check_theory(root):
    root = Path(root).resolve()
    theory = root / 'theory'
    out = theory / 'outputs'
    manifest_path = out / 'theory_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get('schema') == 'sehgal2010_camb_theory_v1' and manifest.get('status') == 'COMPLETE_THEORY_WITH_DECLARED_ASSUMPTIONS', 'Expected a complete theory manifest')
    config = json.loads((theory / 'sehgal2010_config.json').read_text())
    require(manifest['configuration'] == config, 'Theory configuration differs from the packaged input')
    require(manifest['configuration_canonical_sha256'] == canonical_sha256(config),
            'Canonical configuration hash mismatch')
    check_record(theory / 'sehgal2010_config.json', manifest['configuration_file'])
    require(config['expected_camb_version'] == '1.6.1' and manifest['versions']['camb'] == '1.6.1',
            'This audited runtime requires CAMB 1.6.1')
    require(config['published_parameters']['h'] == .71, 'Published Sehgal h must be .71')
    published, completion = config['published_parameters'], config['completion_assumptions']
    require([published[x] for x in ['Omega_b', 'Omega_m', 'ns', 'sigma8_linear_z0']] == [.044, .264, .96, .8],
            'Published parameter vector changed')
    require(completion['r'] == 0 and completion['mnu_eV'] == 0 and completion['num_massive_neutrinos'] == 0,
            'Scalar, massless-neutrino completion changed')
    require(completion['Omega_K'] == 0 and completion['w'] == -1 and completion['wa'] == 0,
            'Flat LambdaCDM completion changed')
    require(config['calculation']['halofit_version'] == 'original', 'Smith2003 completion changed')
    require(manifest['producer_cutoffs'] == config['producer_cutoffs'], 'Cutoff metadata differs')
    cutoffs = config['producer_cutoffs']
    require((cutoffs['nside'], cutoffs['primary_lmax'], cutoffs['kappa_deflection_lmax']) == (2048, 5999, 6143),
            'Native generator cutoff contract changed')
    records = {}
    for name, record in manifest['artifact_manifest'].items():
        require(Path(name).name == name, 'Artifact name is not a basename')
        records[name] = check_record(out / name, record)
    npz_path = out / 'sehgal2010_theory.npz'
    with np.load(npz_path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    ell = arrays['ell']
    require(ell.dtype == np.int64 and np.array_equal(ell, np.arange(6501)), 'Theory multipoles are not zero-based int64 0..6500')
    factor = ell * (ell + 1.) / (2 * np.pi)
    for name in ['unlensed_raw_cl', 'lensed_raw_cl', 'unlensed_D_ell', 'lensed_D_ell', 'producer_unlensed_raw_cl']:
        values = arrays[name]
        require(values.shape == (6501, 4) and values.dtype == np.float64 and np.isfinite(values).all(), f'Invalid {name}')
        contract = manifest['array_contract'][name]
        require(contract['columns'] == ORDER and contract['units'] == 'uK_CMB^2' and contract['shape'] == [6501, 4], f'Incorrect units/order: {name}')
    for kind in ['unlensed', 'lensed']:
        positive_semidefinite(arrays[f'{kind}_raw_cl'], kind)
        np.testing.assert_array_equal(arrays[f'{kind}_D_ell'], arrays[f'{kind}_raw_cl'] * factor[:, None])
    require(np.count_nonzero(arrays['unlensed_raw_cl'][:, 2]) == 0, 'Unlensed scalar BB is not exactly zero')
    primary_cut = cutoffs['primary_lmax']
    np.testing.assert_array_equal(arrays['producer_unlensed_raw_cl'][:primary_cut + 1], arrays['unlensed_raw_cl'][:primary_cut + 1])
    require(np.count_nonzero(arrays['producer_unlensed_raw_cl'][primary_cut + 1:]) == 0, 'Producer primary tail is not exactly zero')
    for name in ['phi_raw_cl', 'phi_PP']:
        values = arrays[name]
        require(values.shape == (6501,) and values.dtype == np.float64 and np.isfinite(values).all() and np.all(values >= 0), f'Invalid {name}')
        require(manifest['array_contract'][name]['units'] == 'dimensionless', f'Incorrect phi units: {name}')
    np.testing.assert_array_equal(arrays['phi_PP'], arrays['phi_raw_cl'] * (ell * (ell + 1.))**2 / (2 * np.pi))
    control_ell = arrays['producer_cutoff_ell']
    require(np.array_equal(control_ell, np.arange(1025)), 'Producer-cutoff output does not cover measurement range')
    for name in ['producer_cutoff_lensed_raw_cl', 'producer_cutoff_lensed_D_ell',
                 'producer_cutoff_sampling2_raw_cl', 'producer_cutoff_sampling2_D_ell', 'producer_cutoff_zero_phi_D_ell']:
        require(arrays[name].shape == (1025, 4) and arrays[name].dtype == np.float64 and np.isfinite(arrays[name]).all(), f'Invalid {name}')
    for prefix in ['producer_cutoff_lensed', 'producer_cutoff_sampling2']:
        np.testing.assert_allclose(arrays[prefix + '_D_ell'], arrays[prefix + '_raw_cl'] * factor[:1025, None], rtol=1e-14, atol=0)
    np.testing.assert_array_equal(arrays['producer_cutoff_zero_phi_D_ell'], arrays['producer_unlensed_raw_cl'][:1025] * factor[:1025, None])
    expected_pp = arrays['phi_PP'].copy()
    expected_pp[6144:] = 0
    np.testing.assert_array_equal(arrays['producer_cutoff_phi_PP'], expected_pp)
    normalization = manifest['normalization']
    expected_As = normalization['trial_As'] * (.8 / normalization['trial_sigma8'])**2
    require(np.isclose(normalization['As_calibrated'], expected_As, rtol=1e-13, atol=0), 'As is not sigma8 normalized by the declared rule')
    residual = normalization['achieved_sigma8_linear_z0'] / .8 - 1
    require(np.isclose(residual, normalization['relative_error'], rtol=0, atol=1e-15), 'Sigma8 residual metadata mismatch')
    require(abs(residual) <= config['calculation']['sigma8_relative_tolerance'], 'Sigma8 normalization exceeds the predeclared tolerance')
    require(normalization['full_solve_count'] == 2 and normalization['same_numerical_settings_in_both_solves'] and normalization['nonlinear_lensing_recalculated_at_final_As'], 'Normalization solve provenance is incomplete')
    source_manifest = theory / config['source_manifest']
    require(sha256(source_manifest) == manifest['source_manifest_sha256'], 'Primary-source manifest hash mismatch')
    source_records = json.loads(source_manifest.read_text())
    require(source_records == manifest['source_records'], 'Source record mismatch')
    source_checks = [check_record(source_manifest.parent / item['filename'], item) for item in source_records]
    comparison = {}
    main = arrays['producer_cutoff_lensed_D_ell']
    alternate = arrays['producer_cutoff_sampling2_D_ell']
    scales = np.abs(main).copy()
    scales[:, 3] = np.sqrt(np.maximum(main[:, 0] * main[:, 1], 0))
    for index, label in enumerate(ORDER):
        valid = (control_ell >= 25) & (control_ell < 1000) & (scales[:, index] > 0)
        comparison[label] = float(np.max(np.abs(alternate[valid, index] - main[valid, index]) / scales[valid, index]))
    return {'status': 'PASS', 'scope': 'Independent input/schema/conversion/hash checks; no CAMB solve or harmonics',
            'theory_manifest_sha256': sha256(manifest_path), 'theory_archive': records['sehgal2010_theory.npz'],
            'checks': {'raw_Cl_units_columns': True, 'D_ell_once': True, 'unlensed_BB_r0': True,
                       'covariance_positive_semidefinite': True, 'phi_PP_normalization': True,
                       'primary5999_phi6143_cutoffs': True, 'zero_phi_identity': True,
                       'sigma8_normalization': True, 'all_input_artifact_hashes': True},
            'derived_input_ombh2': published['Omega_b'] * published['h']**2,
            'derived_input_omch2': (published['Omega_m'] - published['Omega_b']) * published['h']**2,
            'sigma8_relative_error': residual, 'quadrature_1p4_vs2_fractional_auto_or_TE_cross_scale_25_999': comparison,
            'primary_source_checks': source_checks,
            'completion_note': 'Declared tau/neutrino/helium/TCMB/recombination choices define this reference and are not recovered original author inputs.'}


def check_producer_input_spectra(path, metadata, theory_archive):
    """Check transfer algebra using already-produced alms' small spectra."""
    path = Path(path)
    record = check_record(path, metadata['input_spectra'])
    require(metadata['contract_sha256'] == canonical_sha256(metadata['contract']), 'Producer contract hash mismatch')
    master, rid = metadata['master_seed'], metadata['realization_id']
    require(isinstance(master, int) and not isinstance(master, bool) and 0 <= master < 2**64,
            'Invalid uint64 master seed')
    require(isinstance(rid, int) and not isinstance(rid, bool) and 0 <= rid <= 99999, 'Invalid realization ID')
    domain = hashlib.sha256(b'sehgal2010-primary-mmdl-kappa-native2048-v1').digest()
    entropy = np.frombuffer(domain, dtype='<u4').astype(np.uint64).tolist()
    entropy.extend([master % 2**32, master // 2**32, rid % 2**32, rid // 2**32])
    expected_seed = int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint32)[0])
    require(metadata['seed'] == expected_seed, 'Producer task seed does not follow the frozen domain/master/ID policy')
    contract = metadata['contract']
    for name in ['realization_id', 'master_seed', 'seed', 'seed_policy', 'parameters', 'coordinate_provenance']:
        require(contract[name] == metadata[name], f'Producer contract/top-level mismatch: {name}')
    for contract_key, top_key in [('source_kappa_sha256', 'source_kappa'),
                                  ('theory_archive_sha256', 'theory_archive'),
                                  ('theory_manifest_sha256', 'theory_manifest')]:
        require(contract[contract_key] == metadata[top_key]['sha256'], f'Producer contract/source hash mismatch: {contract_key}')
    versions = metadata['versions']
    require(versions['camb'] == '1.6.1' and versions['lenspyx'] == '2.0.0' and versions['ducc0'] in ('0.38', '0.38.0'),
            'Producer remapping/reference runtime differs from the pinned control')
    require(metadata['theory_archive']['sha256'] == sha256(theory_archive), 'Producer theory hash mismatch')
    with np.load(theory_archive, allow_pickle=False) as archive:
        expected_primary = archive['producer_unlensed_raw_cl'][:6000].T.copy()
    with np.load(path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    np.testing.assert_array_equal(arrays['ell_primary'], np.arange(6000))
    np.testing.assert_array_equal(arrays['ell_kappa_phi_deflection'], np.arange(6144))
    np.testing.assert_array_equal(arrays['primary_columns'], ORDER)
    require(int(arrays['realization_id']) == metadata['realization_id'] and int(arrays['seed']) == metadata['seed'], 'Producer input-spectrum ID/seed mismatch')
    require(str(arrays['source_kappa_sha256']) == metadata['source_kappa']['sha256'], 'Saved kappa-source identity mismatch')
    require(str(arrays['theory_archive_sha256']) == metadata['theory_archive']['sha256'], 'Saved theory archive identity mismatch')
    require(str(arrays['theory_manifest_sha256']) == metadata['theory_manifest']['sha256'], 'Saved theory manifest identity mismatch')
    expected_parameters = {'nside': 2048, 'primary_lmax': 5999, 'primary_mmax': 5999,
                           'primary_columns': ORDER, 'primary_units': 'uK_CMB^2', 'synalm_new': True,
                           'kappa_lmax': 6143, 'kappa_mmax': 6143, 'kappa_map2alm_iter': 3,
                           'kappa_map2alm_pol': False, 'kappa_use_weights': False, 'kappa_use_pixel_weights': False,
                           'deflection_lmax': 6143, 'epsilon': 1e-6, 'remapping_pol': True,
                           'beam_applied': False, 'mask_applied': False, 'pixel_window_applied': False,
                           'downgrade_performed': False, 'coordinate_rotation_performed': False}
    for name, value in expected_parameters.items():
        require(metadata['parameters'].get(name) == value, f'Producer parameter changed: {name}')
    require(metadata['theory_phi_used_for_remapping'] is False, 'Continuum CAMB phi must not replace the actual kappa input')
    coordinate = metadata['coordinate_provenance']
    require(coordinate['coordinate_rotation_performed'] is False, 'Unexpected producer coordinate rotation')
    if not coordinate['source_frame_available_in_header']:
        require('owner-declared' in coordinate['selection'] and 'not physically re-certified' in coordinate['caveat'], 'Headerless source-frame uncertainty is not preserved')
    for key in ['unlensed_input_C_ell', 'unlensed_sampled_C_ell']:
        description = metadata['input_spectrum_contract'][key]
        require(description['shape'] == [4, 6000] and description['columns'] == ORDER and description['units'] == 'uK_CMB^2 raw C_l', f'Producer input spectrum units/columns mismatch: {key}')
    for key in ['unlensed_input_C_ell', 'unlensed_sampled_C_ell']:
        require(arrays[key].shape == (4, 6000) and arrays[key].dtype == np.float64, f'Invalid {key}')
        positive_semidefinite(arrays[key].T, key)
        require(np.count_nonzero(arrays[key][2]) == 0, f'Scalar-only sampled/input unlensed BB must be exactly zero: {key}')
    np.testing.assert_array_equal(arrays['unlensed_input_C_ell'], expected_primary)
    for key in ['kappa_C_ell', 'phi_C_ell', 'deflection_C_ell', 'kappa_phi_C_ell',
                'kappa_deflection_C_ell', 'phi_deflection_C_ell', 'phi_transfer', 'deflection_transfer']:
        require(arrays[key].shape == (6144,) and arrays[key].dtype == np.float64 and np.isfinite(arrays[key]).all(), f'Invalid {key}')
    ell = np.arange(6144, dtype=np.float64)
    laplacian = ell * (ell + 1)
    transfer_phi = np.divide(2., laplacian, out=np.zeros_like(ell), where=ell > 0)
    transfer_d = np.divide(2., np.sqrt(laplacian), out=np.zeros_like(ell), where=ell > 0)
    np.testing.assert_array_equal(arrays['phi_transfer'], transfer_phi)
    np.testing.assert_allclose(arrays['deflection_transfer'], transfer_d, rtol=8 * np.finfo(float).eps, atol=0)
    kappa, phi, d = arrays['kappa_C_ell'], arrays['phi_C_ell'], arrays['deflection_C_ell']
    require(np.all(kappa >= 0) and np.all(phi >= 0) and np.all(d >= 0), 'Negative input auto-spectrum')
    # No absolute tolerance can conceal a normalization/sign error at small C_L.
    checks = {'phi_C_ell': kappa * transfer_phi**2,
              'deflection_C_ell': kappa * transfer_d**2,
              'kappa_phi_C_ell': kappa * transfer_phi,
              'kappa_deflection_C_ell': kappa * transfer_d,
              'phi_deflection_C_ell': kappa * transfer_phi * transfer_d}
    for key, expected in checks.items():
        np.testing.assert_allclose(arrays[key], expected, rtol=2e-12, atol=0, err_msg=key)
        require(arrays[key][0] == 0, f'Nonzero monopole in {key}')
    require(transfer_phi[1] == 1 and transfer_d[1] == 2 / np.sqrt(2), 'Dipole handling changed')
    pp = phi * laplacian**2 / (2 * np.pi)
    require(arrays['phi_PP'].shape == (6144,) and arrays['phi_PP'].dtype == np.float64 and np.isfinite(arrays['phi_PP']).all(), 'Invalid saved phi_PP')
    np.testing.assert_array_equal(arrays['phi_PP'], pp)
    np.testing.assert_allclose(pp[1:], 2 * kappa[1:] / np.pi, rtol=2e-12, atol=0)
    require(pp[0] == 0, 'Phi PP monopole must vanish')
    report = {'status': 'PASS', 'realization_id': metadata['realization_id'], 'input_spectra': record,
              'theory_sha256': metadata['theory_archive']['sha256'],
              'checks': {'unlensed_raw_Cl_identity_and_order': True, 'sampled_unlensed_BB_exact_zero': True,
                         'kappa_phi_deflection_auto_cross_sign_normalization': True,
                         'phi_monopole_zero_dipole_retained': True, 'dimensionless_PP_conversion': True},
              'scope': 'Checks existing producer algebra/inputs; no independent kappa-map harmonic transform or physical-frame certification'}
    return report, arrays, pp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = check_theory(args.root)
    rendered = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
