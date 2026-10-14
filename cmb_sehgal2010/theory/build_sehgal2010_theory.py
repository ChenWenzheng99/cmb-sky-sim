#!/usr/bin/env python3
"""Build explicit-completion Sehgal2010 CAMB primary and lensed spectra.

This standalone builder imports only CAMB, NumPy and the standard library.
It does not read or overwrite any old CMB map, FFP10 input or release file.
"""
import argparse
import gc
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import platform
import tempfile
import time


SCHEMA = 'sehgal2010_camb_theory_v1'
COLUMNS = ['TT', 'EE', 'BB', 'TE']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def quadrature_comparison(np, main_D, check_D):
    comparison = {}
    for i, name in enumerate(COLUMNS):
        main = main_D[2:, i]
        check = check_D[2:, i]
        nonzero = main != 0
        comparison[name] = {
            'ell_min': 2, 'ell_max': len(main_D) - 1,
            'max_abs_D_ell_difference': float(np.max(np.abs(check - main))),
            'max_abs_relative_difference': float(np.max(np.abs(check[nonzero] / main[nonzero] - 1)))
            if np.any(nonzero) else None}
    return comparison


def canonicalize_correlation_spectra(np, returned, ell):
    """Retain CAMB's special zero mode, export canonical raw C/D with ell0=0.

    CAMB1.6.1 correlations internally sets lfacs[0]=1. Its returned ell0
    therefore is not physical ell(ell+1)C_ell/(2pi). Positive ell values
    are left unchanged; the nonstandard zero mode is retained separately.
    """
    if returned.shape != (len(ell), 4) or not np.isfinite(returned).all():
        raise ValueError('Invalid correlations return shape/finite domain')
    raw_zero_mode = returned[0].copy()
    canonical_D = returned.copy()
    canonical_D[0] = 0.0
    factor = ell * (ell + 1) / (2 * np.pi)
    raw = np.divide(canonical_D, factor[:, None], out=np.zeros_like(canonical_D),
                    where=factor[:, None] > 0)
    return raw, canonical_D, raw_zero_mode


def load_config(path):
    path = Path(path).resolve()
    configuration = json.loads(path.read_text())
    if configuration.get('schema') != SCHEMA:
        raise ValueError(f'Expected schema {SCHEMA}')
    published = configuration['published_parameters']
    assumed = configuration['completion_assumptions']
    calc = configuration['calculation']
    producer = configuration['producer_cutoffs']
    for name in ['Omega_b', 'Omega_m', 'h', 'ns', 'sigma8_linear_z0']:
        value = published[name]
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'Invalid {name}')
    if published['Omega_m'] <= published['Omega_b']:
        raise ValueError('Omega_m must exceed Omega_b')
    for name in ['table_lmax', 'lens_margin', 'lens_potential_accuracy']:
        value = calc[name]
        if isinstance(value, bool) or not isinstance(value, int) or value < (2 if name == 'table_lmax' else 0):
            raise ValueError(f'Invalid {name}')
    for name in ['accuracy_boost', 'l_accuracy_boost', 'l_sample_boost',
                 'k_eta_max_per_internal_ell', 'matter_kmax_Mpc_inverse',
                 'trial_As', 'sigma8_relative_tolerance']:
        value = calc[name]
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'Invalid {name}')
    for name in ['primary_lmax', 'kappa_deflection_lmax', 'nside']:
        value = producer[name]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise ValueError(f'Invalid producer {name}')
    if producer['kappa_deflection_lmax'] != 3 * producer['nside'] - 1:
        raise ValueError('Deflection cutoff must match the retained native-grid contract')
    if calc['table_lmax'] < max(producer['primary_lmax'], producer['kappa_deflection_lmax']):
        raise ValueError('Theoretical table must cover both producer cutoffs')
    if assumed['r'] != 0 or assumed['mnu_eV'] != 0 or assumed['num_massive_neutrinos'] != 0:
        raise ValueError('This release defines scalar-only, massless-neutrino completion')
    if assumed['Omega_K'] != 0 or assumed['w'] != -1 or assumed['wa'] != 0:
        raise ValueError('This release defines flat LambdaCDM')
    if calc['halofit_version'] != 'original':
        raise ValueError('This release explicitly uses original Smith2003 nonlinear lensing')
    control_lmax = calc['cutoff_control_output_lmax']
    if isinstance(control_lmax, bool) or not isinstance(control_lmax, int) or not 2 <= control_lmax <= calc['table_lmax']:
        raise ValueError('Invalid cutoff control output lmax')
    control_sampling = calc['cutoff_control_sampling_factor']
    if isinstance(control_sampling, bool) or not isinstance(control_sampling, (float, int)) or not math.isfinite(control_sampling) or control_sampling <= 0:
        raise ValueError('Invalid cutoff control sampling factor')
    check_sampling = calc['cutoff_control_check_sampling_factor']
    if isinstance(check_sampling, bool) or not isinstance(check_sampling, (float, int)) or not math.isfinite(check_sampling) or check_sampling <= control_sampling:
        raise ValueError('Cutoff control check sampling factor must exceed the main value')
    if not isinstance(calc['accurate_BB'], bool):
        raise ValueError('accurate_BB must be boolean')
    return path, configuration


def make_parameters(camb, configuration, scalar_amplitude):
    published = configuration['published_parameters']
    assumed = configuration['completion_assumptions']
    calc = configuration['calculation']
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=100 * published['h'],
                      ombh2=published['Omega_b'] * published['h'] ** 2,
                      omch2=(published['Omega_m'] - published['Omega_b']) * published['h'] ** 2,
                      omk=assumed['Omega_K'], tau=assumed['tau'],
                      mnu=assumed['mnu_eV'], num_massive_neutrinos=assumed['num_massive_neutrinos'],
                      nnu=assumed['N_eff'], YHe=assumed['YHe'], TCMB=assumed['TCMB_K'])
    pars.set_dark_energy(w=assumed['w'], wa=assumed['wa'], dark_energy_model='fluid')
    pars.InitPower.set_params(As=scalar_amplitude, ns=published['ns'], r=assumed['r'],
                              nrun=assumed['nrun'], nrunrun=assumed['nrunrun'],
                              pivot_scalar=assumed['pivot_scalar_Mpc_inverse'])
    pars.WantTensors = False
    pars.set_accuracy(AccuracyBoost=calc['accuracy_boost'],
                      lAccuracyBoost=calc['l_accuracy_boost'], lSampleBoost=calc['l_sample_boost'])
    pars.Accuracy.AccurateBB = calc['accurate_BB']
    pars.set_matter_power(redshifts=[0.0], kmax=calc['matter_kmax_Mpc_inverse'])
    internal_lmax = calc['table_lmax'] + calc['lens_margin']
    pars.set_for_lmax(calc['table_lmax'], lens_margin=calc['lens_margin'],
                      lens_potential_accuracy=calc['lens_potential_accuracy'],
                      max_eta_k=calc['k_eta_max_per_internal_ell'] * internal_lmax)
    pars.NonLinear = camb.model.NonLinear_both
    pars.NonLinearModel.set_params(halofit_version=calc['halofit_version'])
    pars.Alens = 1.0
    return pars


def validate_arrays(np, arrays, configuration):
    limit = configuration['calculation']['table_lmax']
    ell = arrays['ell']
    np.testing.assert_array_equal(ell, np.arange(limit + 1))
    if ell.dtype != np.int64:
        raise ValueError('ell must be int64')
    for key in ['unlensed_raw_cl', 'lensed_raw_cl', 'unlensed_D_ell', 'lensed_D_ell',
                'producer_unlensed_raw_cl']:
        if arrays[key].shape != (limit + 1, 4) or arrays[key].dtype != np.float64 or not np.isfinite(arrays[key]).all():
            raise ValueError(f'Invalid {key} shape/dtype/finite domain')
    if np.count_nonzero(arrays['unlensed_raw_cl'][:, 2]):
        raise ValueError('Scalar-only unlensed BB must be exactly zero')
    for key in ['unlensed_raw_cl', 'lensed_raw_cl']:
        values = arrays[key]
        if np.any(values[:, :3] < 0):
            raise ValueError(f'Negative auto-spectrum in {key}')
        violation = values[:, 3] ** 2 - values[:, 0] * values[:, 1]
        if np.any(violation > 1e-10 * np.maximum(values[:, 0] * values[:, 1], 1e-300)):
            raise ValueError(f'T/E covariance is not positive semidefinite in {key}')
    factor = ell * (ell + 1) / (2 * np.pi)
    np.testing.assert_array_equal(arrays['unlensed_D_ell'], arrays['unlensed_raw_cl'] * factor[:, None])
    np.testing.assert_array_equal(arrays['lensed_D_ell'], arrays['lensed_raw_cl'] * factor[:, None])
    for key in ['phi_raw_cl', 'phi_PP']:
        if arrays[key].shape != (limit + 1,) or not np.isfinite(arrays[key]).all() or np.any(arrays[key] < 0):
            raise ValueError(f'Invalid {key}')
    np.testing.assert_array_equal(arrays['phi_PP'], arrays['phi_raw_cl'] * (ell * (ell + 1)) ** 2 / (2 * np.pi))
    cutoff = configuration['producer_cutoffs']['primary_lmax']
    np.testing.assert_array_equal(arrays['producer_unlensed_raw_cl'][:cutoff + 1], arrays['unlensed_raw_cl'][:cutoff + 1])
    if np.count_nonzero(arrays['producer_unlensed_raw_cl'][cutoff + 1:]):
        raise ValueError('Producer primary spectrum must vanish above its cutoff')
    control_limit = configuration['calculation']['cutoff_control_output_lmax']
    np.testing.assert_array_equal(arrays['producer_cutoff_ell'], np.arange(control_limit + 1, dtype=np.int64))
    for key in ['producer_cutoff_lensed_raw_cl', 'producer_cutoff_lensed_D_ell',
                'producer_cutoff_sampling2_raw_cl', 'producer_cutoff_sampling2_D_ell',
                'producer_cutoff_zero_phi_D_ell']:
        if arrays[key].shape != (control_limit + 1, 4) or arrays[key].dtype != np.float64 or not np.isfinite(arrays[key]).all():
            raise ValueError(f'Invalid {key}')
    np.testing.assert_allclose(arrays['producer_cutoff_lensed_D_ell'],
                               arrays['producer_cutoff_lensed_raw_cl'] * factor[:control_limit + 1, None],
                               rtol=1e-14, atol=0)
    expected_cutoff_pp = arrays['phi_PP'].copy()
    expected_cutoff_pp[configuration['producer_cutoffs']['kappa_deflection_lmax'] + 1:] = 0.0
    np.testing.assert_array_equal(arrays['producer_cutoff_phi_PP'], expected_cutoff_pp)
    if arrays['raw_correlation_zero_mode'].shape != (3, 4) or not np.isfinite(arrays['raw_correlation_zero_mode']).all():
        raise ValueError('Invalid retained correlation zero modes')
    for key in ['producer_cutoff_lensed_raw_cl', 'producer_cutoff_lensed_D_ell',
                'producer_cutoff_sampling2_raw_cl', 'producer_cutoff_sampling2_D_ell',
                'producer_cutoff_zero_phi_D_ell']:
        if np.count_nonzero(arrays[key][0]):
            raise ValueError(f'Canonical ell0 must vanish: {key}')
    np.testing.assert_array_equal(arrays['producer_cutoff_zero_phi_D_ell'],
                                   arrays['producer_unlensed_raw_cl'][:control_limit + 1] * factor[:control_limit + 1, None])
    return {'ell_index_exact': True, 'float64_finite': True, 'unlensed_BB_exact_zero': True,
            'auto_spectra_nonnegative': True, 'TE_covariance_positive_semidefinite': True,
            'C_ell_D_ell_conversion_exact': True, 'phi_PP_conversion_exact': True,
            'producer_cutoff_exact': True, 'cutoff_control_finite': True,
            'cutoff_control_C_ell_D_ell_roundtrip': True, 'cutoff_control_phi_cutoff_exact': True,
            'zero_phi_D_ell_identity_exact': True, 'canonical_auxiliary_ell0_zero': True,
            'nonstandard_correlations_zero_modes_retained': True}


def build(config_path, configuration, out, overwrite=False):
    import camb
    from camb import correlations
    import numpy as np

    if camb.__version__ != configuration['expected_camb_version']:
        raise ValueError(f"CAMB version {camb.__version__} differs from frozen {configuration['expected_camb_version']}")
    out = Path(out).resolve()
    output_names = ['sehgal2010_theory.npz', 'sehgal2010_unlensed_rawCl.dat',
                    'sehgal2010_lensed_rawCl.dat', 'theory_manifest.json']
    if not overwrite and any((out / name).exists() for name in output_names):
        raise FileExistsError(f'Theory output already exists in {out}; choose a fresh output directory')
    source_manifest_path = config_path.parent / configuration['source_manifest']
    source_records = json.loads(source_manifest_path.read_text())
    for item in source_records:
        if sha256(source_manifest_path.parent / item['filename']) != item['sha256']:
            raise ValueError(f"Primary-source snapshot hash mismatch: {item['filename']}")
    started = time.monotonic()
    target = configuration['published_parameters']['sigma8_linear_z0']
    amplitude_trial = configuration['calculation']['trial_As']
    trial_pars = make_parameters(camb, configuration, amplitude_trial)
    print('CAMB full solve 1/2: fixed-precision sigma8 normalization trial', flush=True)
    trial = camb.get_results(trial_pars)
    sigma_trial = float(trial.get_sigma8_0())
    if not math.isfinite(sigma_trial) or sigma_trial <= 0:
        raise ValueError('Invalid trial sigma8')
    amplitude_final = amplitude_trial * (target / sigma_trial) ** 2
    del trial
    gc.collect()
    pars = make_parameters(camb, configuration, amplitude_final)
    print('CAMB full solve 2/2: final As with nonlinear lensing recomputed', flush=True)
    results = camb.get_results(pars)
    sigma_final = float(results.get_sigma8_0())
    relative_error = sigma_final / target - 1
    if abs(relative_error) > configuration['calculation']['sigma8_relative_tolerance']:
        raise ValueError(f'Sigma8 normalization failed: target={target}, achieved={sigma_final}')
    lmax = configuration['calculation']['table_lmax']
    ell = np.arange(lmax + 1, dtype=np.int64)
    unlensed = results.get_unlensed_scalar_cls(lmax=lmax, CMB_unit='muK', raw_cl=True)
    lensed = results.get_lensed_scalar_cls(lmax=lmax, CMB_unit='muK', raw_cl=True)
    phi = results.get_lens_potential_cls(lmax=lmax, raw_cl=True)[:, 0]
    producer = unlensed.copy()
    producer[configuration['producer_cutoffs']['primary_lmax'] + 1:] = 0.0
    factor = ell * (ell + 1) / (2 * np.pi)
    arrays = {'ell': ell, 'unlensed_raw_cl': unlensed, 'lensed_raw_cl': lensed,
              'unlensed_D_ell': unlensed * factor[:, None], 'lensed_D_ell': lensed * factor[:, None],
              'phi_raw_cl': phi, 'phi_PP': phi * (ell * (ell + 1)) ** 2 / (2 * np.pi),
              'producer_unlensed_raw_cl': producer}
    control_limit = configuration['calculation']['cutoff_control_output_lmax']
    cutoff_pp = arrays['phi_PP'].copy()
    cutoff_pp[configuration['producer_cutoffs']['kappa_deflection_lmax'] + 1:] = 0.0
    print('Full-theta producer-cutoff reference: sampling factor 1.4', flush=True)
    control_D = correlations.lensed_cls(producer * factor[:, None], cutoff_pp,
                                       lmax=lmax, lmax_lensed=control_limit, theta_max=None,
                                       sampling_factor=configuration['calculation']['cutoff_control_sampling_factor'],
                                       leggaus=True, cache=True)
    control_raw, control_D, control_zero_mode = canonicalize_correlation_spectra(np, control_D, ell[:control_limit + 1])
    print('Full-theta quadrature comparison: sampling factor 2.0', flush=True)
    sampling2_D = correlations.lensed_cls(producer * factor[:, None], cutoff_pp,
                                         lmax=lmax, lmax_lensed=control_limit, theta_max=None,
                                         sampling_factor=configuration['calculation']['cutoff_control_check_sampling_factor'],
                                         leggaus=True, cache=True)
    sampling2_raw, sampling2_D, sampling2_zero_mode = canonicalize_correlation_spectra(np, sampling2_D, ell[:control_limit + 1])
    print('Full-theta zero-phi exact-identity control', flush=True)
    zero_phi_D = correlations.lensed_cls(producer * factor[:, None], np.zeros_like(cutoff_pp),
                                        lmax=lmax, lmax_lensed=control_limit, theta_max=None,
                                        sampling_factor=configuration['calculation']['cutoff_control_sampling_factor'],
                                        leggaus=True, cache=True)
    _, zero_phi_D, zero_phi_zero_mode = canonicalize_correlation_spectra(np, zero_phi_D, ell[:control_limit + 1])
    arrays.update(producer_cutoff_ell=np.arange(control_limit + 1, dtype=np.int64),
                  producer_cutoff_phi_PP=cutoff_pp,
                  producer_cutoff_lensed_D_ell=control_D,
                  producer_cutoff_lensed_raw_cl=control_raw,
                  producer_cutoff_sampling2_raw_cl=sampling2_raw,
                  producer_cutoff_sampling2_D_ell=sampling2_D,
                  producer_cutoff_zero_phi_D_ell=zero_phi_D,
                  raw_correlation_zero_mode=np.stack([control_zero_mode, sampling2_zero_mode, zero_phi_zero_mode]))
    checks = validate_arrays(np, arrays, configuration)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.sehgal2010_theory_', dir=out) as stage_name:
        stage = Path(stage_name)
        np.savez_compressed(stage / output_names[0], **arrays)
        for name, values in zip(output_names[1:3], [unlensed, lensed]):
            np.savetxt(stage / name, np.column_stack([ell, values]), fmt=['%d'] + ['%.17e'] * 4,
                       header='ell C_TT C_EE C_BB C_TE; raw C_ell, thermodynamic uK_CMB^2; zero-based ell')
        array_contract = {key: {'shape': list(value.shape), 'dtype': str(value.dtype),
                               'units': ('nonstandard_correlation_zero_mode_uK2' if key == 'raw_correlation_zero_mode' else
                                         'dimensionless' if key in ['phi_raw_cl', 'phi_PP', 'ell', 'producer_cutoff_ell', 'producer_cutoff_phi_PP'] else 'uK_CMB^2'),
                               **({'columns': COLUMNS} if value.ndim == 2 else {})}
                          for key, value in arrays.items()}
        manifest = {
            'schema': SCHEMA, 'status': 'COMPLETE_THEORY_WITH_DECLARED_ASSUMPTIONS',
            'configuration': configuration, 'configuration_canonical_sha256': canonical_sha256(configuration),
            'configuration_file': {'path': str(config_path), 'sha256': sha256(config_path)},
            'producer_cutoffs': configuration['producer_cutoffs'], 'array_contract': array_contract,
            'correlation_zero_mode_convention': {
                'canonical_auxiliary_raw_C0_and_D0': 0.0,
                'source_behavior': 'CAMB1.6.1 correlations uses lfacs[0]=1; returned ell0 is not physical D_ell',
                'raw_correlation_zero_mode': arrays['raw_correlation_zero_mode'].tolist(),
                'rows': ['producer_cutoff_sampling1.4', 'producer_cutoff_sampling2.0', 'zero_phi_sampling1.4'],
                'columns': COLUMNS, 'positive_ell_correlation_output_unchanged': True,
                'main_CAMB_arrays_unchanged': True},
            'producer_cutoff_control': {'primary_lmax': configuration['producer_cutoffs']['primary_lmax'],
                                        'phi_lmax': configuration['producer_cutoffs']['kappa_deflection_lmax'],
                                        'output_lmax': control_limit, 'angular_range': 'full [0, pi]',
                                        'sampling_factor': configuration['calculation']['cutoff_control_sampling_factor'],
                                        'check_sampling_factor': configuration['calculation']['cutoff_control_check_sampling_factor'],
                                        'quadrature_comparison': quadrature_comparison(np, control_D, sampling2_D),
                                        'phi_source': 'Self-consistent CAMB phi_PP, exact zero above producer cutoff',
                                        'purpose': 'Bandlimit-matched continuum control; the actual mmDL kappa map is not substituted by CAMB phi'},
            'normalization': {'target_sigma8_linear_z0': target, 'trial_As': amplitude_trial,
                              'trial_sigma8': sigma_trial, 'As_calibrated': amplitude_final,
                              'achieved_sigma8_linear_z0': sigma_final, 'relative_error': relative_error,
                              'full_solve_count': 2, 'same_numerical_settings_in_both_solves': True,
                              'nonlinear_lensing_recalculated_at_final_As': True},
            'camb_actual_parameters_string': str(results.Params),
            'derived_parameters': {key: float(value) for key, value in results.get_derived_params().items()},
            'background_Omega_z0': {name: float(results.get_Omega(name))
                                   for name in ['K', 'cdm', 'baryon', 'photon', 'neutrino', 'nu', 'de']},
            'checks': checks, 'lmax_internal_actual': int(results.Params.max_l),
            'max_eta_k_actual': float(results.Params.max_eta_k),
            'versions': {'camb': camb.__version__, 'numpy': np.__version__, 'python': platform.python_version()},
            'source_records': source_records, 'source_manifest_sha256': sha256(source_manifest_path),
            'builder_sha256': sha256(__file__),
            'installed_camb_source_sha256': {name: sha256(inspect.getsourcefile(module))
                                            for name, module in [('model', camb.model), ('results', camb.results), ('correlations', correlations)]},
            'artifact_manifest': {name: {'path': str(out / name), 'bytes': (stage / name).stat().st_size,
                                        'sha256': sha256(stage / name)} for name in output_names[:3]},
            'host': platform.node(), 'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'elapsed_seconds': time.monotonic() - started,
            'scope': 'Self-consistent CAMB primary and continuum nonlinear-lensing reference for the stated completion. Producer remaps the scalar primary with separately stored mmDL kappa; CAMB phi is not substituted for that map.',
            'limitations': ['Complete original Sehgal CAMB ini and exact CMB tau are not recovered.',
                            'The scalar-only completion sets primordial BB to zero.',
                            'Full CAMB lensed arrays retain the calculated high-ell power; producer arrays explicitly cut primary at 5999.',
                            'CAMB continuum phi does not include the released map finite-box, replication, GAN or particle effects.',
                            'Numerical settings are declared; a resolution/accuracy convergence claim requires separate controls.']}
        (stage / output_names[3]).write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n')
        for name in output_names:
            os.replace(stage / name, out / name)
    print(f'COMPLETE {out / output_names[0]} sigma8={sigma_final:.12g} As={amplitude_final:.12g}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('sehgal2010_config.json'))
    parser.add_argument('--out', type=Path, help='Fresh theory-output directory')
    parser.add_argument('--print-config', action='store_true', help='Validate/print without importing or running CAMB')
    parser.add_argument('--overwrite', action='store_true', help='Replace files only in this selected theory-output directory')
    args = parser.parse_args()
    path, configuration = load_config(args.config)
    if args.print_config:
        print(json.dumps(configuration, indent=2, allow_nan=False))
        return
    out = args.out or path.parent / configuration['output_dir']
    build(path, configuration, out, args.overwrite)


if __name__ == '__main__':
    main()
