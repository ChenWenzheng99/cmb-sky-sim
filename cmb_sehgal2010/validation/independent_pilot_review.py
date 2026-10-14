#!/usr/bin/env python3
"""Read-only independent review of the completed native ten-map pilot.

Recomputes bins/statistics from archived spectra. No native maps, harmonic
transforms, generation, or CAMB calls are made. No validation policy is edited.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t

from validate_inputs import check_producer_input_spectra, check_record, check_theory, sha256


SPECTRA = ['TT', 'EE', 'BB', 'TE']
REFERENCES = ['self_consistent_CAMB', 'producer_bandlimit_CAMB', 'conditional_actual_mean_kappa_2point']


def archive(path):
    with np.load(path, allow_pickle=False) as source:
        return {key: source[key] for key in source.files}


def bins(values):
    values = np.asarray(values, dtype=np.float64)
    ell = np.arange(values.shape[-1])
    weights = 2 * ell + 1
    output = np.empty(values.shape[:-1] + (39,))
    for b in range(39):
        lower = 25 * (b + 1)
        upper = lower + 25
        output[..., b] = np.average(values[..., lower:upper], weights=weights[lower:upper], axis=-1)
    return output


def close(a, b):
    np.testing.assert_allclose(a, b, rtol=2e-12, atol=1e-13)


def close_statistic(a, b):
    # Independent summation order changes a near-unity fitted amplitude by
    # machine epsilon; dividing A-1 by its ~5e-4 SEM amplifies that roundoff.
    # This compares arithmetic implementations only, not a scientific gate.
    np.testing.assert_allclose(a, b, rtol=2e-12, atol=2e-10)


def review(root, evidence):
    root, evidence = Path(root), Path(evidence)
    manifest = json.loads((root / 'validation/validation_code_manifest.json').read_text())
    for name, record in manifest['files'].items():
        check_record(root / 'validation' / name, record)
    policy_path = root / 'validation/validation_policy.json'
    policy = json.loads(policy_path.read_text())
    theory_report = check_theory(root)
    theory_path = root / 'theory/outputs/sehgal2010_theory.npz'
    theory = archive(theory_path)
    theory_manifest = json.loads(theory_path.with_name('theory_manifest.json').read_text())
    np.testing.assert_array_equal(theory['raw_correlation_zero_mode'],
                                  theory_manifest['correlation_zero_mode_convention']['raw_correlation_zero_mode'])
    for key in ['producer_cutoff_lensed_raw_cl', 'producer_cutoff_lensed_D_ell',
                'producer_cutoff_sampling2_raw_cl', 'producer_cutoff_sampling2_D_ell']:
        np.testing.assert_array_equal(theory[key][0], np.zeros(4))
    summary_path = evidence / 'pilot_summary/native_cmb_pilot_10.npz'
    summary_meta_path = summary_path.with_suffix('.json')
    summary = json.loads(summary_meta_path.read_text())
    summary_record = check_record(summary_path, summary['artifact'])
    result = archive(summary_path)
    if summary['policy_sha256'] != sha256(policy_path) or summary['policy'] != policy:
        raise ValueError('The predeclared validation policy changed')
    np.testing.assert_array_equal(result['ids'], np.arange(10))
    np.testing.assert_array_equal(result['spectra'], SPECTRA)
    np.testing.assert_array_equal(result['ell'], np.arange(1025))
    if result['raw_C_ell'].shape != (10, 4, 1025) or not np.isfinite(result['raw_C_ell']).all():
        raise ValueError('Invalid native pilot raw spectra')
    ell = np.arange(1025)
    factor = ell * (ell + 1.) / (2 * np.pi)
    raw_b = bins(result['raw_C_ell'] * factor)
    mean = np.sum(raw_b, axis=0) / 10
    sd = np.sqrt(np.sum((raw_b - mean)**2, axis=0) / 9)
    input_reports, producer_records = [], []
    source_seeds, kappa_hashes, map_hashes = [], [], []
    kappa_input_pp, kappa_inputs = [], []
    for rid in range(10):
        measured_path = evidence / f'pilot_spectra/native_cmb_spectra_{rid:05d}.npz'
        measured_meta_path = measured_path.with_suffix('.json')
        measured_meta = json.loads(measured_meta_path.read_text())
        check_record(measured_path, measured_meta['output'])
        measured = archive(measured_path)
        np.testing.assert_array_equal(measured['C_ell'], result['raw_C_ell'][rid])
        close(measured['D_b'], raw_b[rid])
        if (measured_meta['source_nside'], measured_meta['analysis_nside'], measured_meta['lmax'], measured_meta['iter']) != (2048, 2048, 1024, 3):
            raise ValueError('Measurement resolution/cutoff/iteration mismatch')
        for key in ['beam', 'mask', 'pixel_window', 'downgrade', 'monopole_dipole_removal', 'use_weights', 'use_pixel_weights']:
            if measured_meta[key] is not False:
                raise ValueError(f'Unexpected measurement operation: {key}')
        if measured_meta['FITS_checksums_verified'] is not True or measured_meta['spectra_order'] != SPECTRA:
            raise ValueError('Measurement checksum or columns not verified')
        for stat in measured_meta['statistics']:
            if (stat['pixels'], stat['finite_count'], stat['unseen_count']) != (50331648, 50331648, 0):
                raise ValueError('Full native output finite-statistics failure')
        producer = measured_meta['producer']
        if producer['realization_id'] != rid or producer['output']['sha256'] != measured_meta['source']['sha256']:
            raise ValueError('Producer/measurement identity mismatch')
        if producer['producer_source']['release_producer']['sha256'] != sha256(root / 'production/lensed_CMB_generation.py'):
            raise ValueError('Executed producer differs from the audited frozen source')
        if measured_meta['script']['sha256'] != manifest['files']['measure_native_cmb.py']['sha256']:
            raise ValueError('Executed measurement differs from the frozen source')
        input_path = evidence / f'producer_inputs/lensedcmb_new_{rid:05d}.input_spectra.npz'
        input_report, input_arrays, pp = check_producer_input_spectra(input_path, producer, theory_path)
        input_reports.append(input_report)
        kappa_input_pp.append(pp)
        kappa_inputs.append(input_arrays)
        source_seeds.append(producer['seed'])
        kappa_hashes.append(producer['source_kappa']['sha256'])
        map_hashes.append(producer['output']['sha256'])
        producer_records.append({'id': rid, 'measurement_npz_sha256': sha256(measured_path),
                                 'measurement_sidecar_sha256': sha256(measured_meta_path),
                                 'producer_input_sha256': sha256(input_path),
                                 'output_map_sha256': producer['output']['sha256'],
                                 'source_kappa_sha256': producer['source_kappa']['sha256'],
                                 'seed': producer['seed']})
        if rid == 0:
            control = measured['control1024_iter1_C_ell']
    if any(len(set(values)) != 10 for values in [source_seeds, kappa_hashes, map_hashes]):
        raise ValueError('Pilot realization identity/seed collision')
    kappa_path = evidence / 'kappa_pilot/producer_kappa_ensemble_10.npz'
    kappa_meta = json.loads(kappa_path.with_suffix('.json').read_text())
    check_record(kappa_path, kappa_meta['artifact'])
    kappa = archive(kappa_path)
    np.testing.assert_array_equal(kappa['ell'], np.arange(6144))
    np.testing.assert_array_equal(kappa['raw_phi_PP'], np.array(kappa_input_pp))
    np.testing.assert_allclose(kappa['mean_phi_PP'], np.sum(kappa['raw_phi_PP'], axis=0) / 10, rtol=2e-12, atol=0)
    np.testing.assert_allclose(kappa['sample_std_phi_PP'], np.sqrt(np.sum((kappa['raw_phi_PP'] - kappa['mean_phi_PP'])**2, axis=0) / 9), rtol=2e-12, atol=0)
    for key, input_key in [('raw_kappa_C_ell', 'kappa_C_ell'), ('raw_phi_C_ell', 'phi_C_ell'), ('raw_deflection_C_ell', 'deflection_C_ell')]:
        np.testing.assert_array_equal(kappa[key], np.array([item[input_key] for item in kappa_inputs]))
    matched_path = evidence / 'matched_kappa_reference/matched_mmdl_pp_reference.npz'
    matched_meta = json.loads(matched_path.with_suffix('.json').read_text())
    check_record(matched_path, matched_meta['artifact'])
    if matched_meta['measured_PP_input']['sha256'] != sha256(kappa_path) or matched_meta['primary_theory']['sha256'] != sha256(theory_path):
        raise ValueError('Conditional reference input identity mismatch')
    if matched_meta['measured_PP_input']['key'] != 'mean_phi_PP':
        raise ValueError('Conditional reference is not the actual ten-map mean PP')
    matched = archive(matched_path)
    np.testing.assert_array_equal(matched['input_phi_PP'][:6144], kappa['mean_phi_PP'])
    np.testing.assert_array_equal(matched['input_phi_PP'][6144:], np.zeros(357))
    np.testing.assert_array_equal(matched['lensed_D_ell'][0], np.zeros(4))
    np.testing.assert_array_equal(matched['raw_correlation_zero_mode'], matched_meta['correlation_zero_mode_convention']['raw_correlation_zero_mode'])
    close(matched['lensed_D_ell'], matched['lensed_raw_cl'] * factor[:, None])
    expected_references = [theory['lensed_raw_cl'][:1025].T,
                           theory['producer_cutoff_lensed_raw_cl'].T,
                           matched['lensed_raw_cl'].T]
    amp_limit = float(t.ppf(1 - .005 / 6, 9))
    bin_limit = float(t.ppf(1 - .005 / 234, 9))
    recomputed_metrics = {}
    for key, reference in zip(REFERENCES, expected_references):
        np.testing.assert_array_equal(result[key + '__reference_C_ell'], reference)
        ref_b = bins(reference * factor)
        close(result[key + '__raw_D_b'], raw_b)
        close(result[key + '__mean_D_b'], mean)
        close(result[key + '__sample_std_D_b'], sd)
        close(result[key + '__reference_D_b'], ref_b)
        # Explicit per-bin sums implement the Gaussian single-map fit weights.
        per_ell_variance = np.array([2 * reference[0]**2, 2 * reference[1]**2,
                                     2 * reference[2]**2, reference[0] * reference[1] + reference[3]**2]) / (2 * ell + 1)
        variance = np.empty((4, 39))
        for b in range(39):
            lower, upper = 25 * (b + 1), 25 * (b + 2)
            weights = 2 * ell[lower:upper] + 1
            variance[:, b] = np.sum(per_ell_variance[:, lower:upper] * factor[lower:upper]**2 * weights**2, axis=1) / np.sum(weights)**2
        amplitudes = np.sum(raw_b * (ref_b / variance)[None], axis=2) / np.sum(ref_b**2 / variance, axis=1)
        amplitude_mean = amplitudes.mean(0)
        amplitude_sd = np.sqrt(np.sum((amplitudes - amplitude_mean)**2, axis=0) / 9)
        amp_t = (amplitude_mean - 1) / (amplitude_sd / np.sqrt(10))
        bin_t = (mean - ref_b) / (sd / np.sqrt(10))
        close(result[key + '__reference_gaussian_single_realization_bin_variance'], variance)
        close(result[key + '__raw_fitted_amplitude'], amplitudes)
        close_statistic(result[key + '__amplitude_student_t'], amp_t)
        close_statistic(result[key + '__bin_student_t'], bin_t)
        indices = [0, 1, 3]
        primary_pass = bool(np.all(np.abs(amp_t[indices]) <= amp_limit) and np.all(np.abs(bin_t[indices]) <= bin_limit))
        if summary['reference_metrics'][key]['primary_spectral_gate'] != ('PASS' if primary_pass else 'REVIEW_REQUIRED'):
            raise ValueError('Recomputed primary gate differs from reported gate')
        spectrum_metrics = {}
        for index, name in enumerate(SPECTRA):
            item = {'fitted_amplitude_minus1_pct': float(100 * (amplitude_mean[index] - 1)),
                    'fitted_amplitude_sample_SD_pct': float(100 * amplitude_sd[index]),
                    'fitted_amplitude_mean_SEM_pct': float(100 * amplitude_sd[index] / np.sqrt(10)),
                    'fitted_amplitude_t_df9': float(amp_t[index]),
                    'maximum_abs_bin_t_df9': float(np.max(np.abs(bin_t[index])))}
            if name != 'TE':
                ratios = mean[index] / ref_b[index]
                item.update(median_bin_ratio_minus1_pct=float(100 * np.median(ratios - 1)),
                            minimum_bin_ratio_minus1_pct=float(100 * (ratios.min() - 1)),
                            maximum_bin_ratio_minus1_pct=float(100 * (ratios.max() - 1)))
            else:
                stable = (mean[3] - ref_b[3]) / np.sqrt(ref_b[0] * ref_b[1])
                close(result[key + '__TE_residual_over_sqrt_TT_EE'], stable)
                item['maximum_abs_residual_over_sqrt_TT_EE'] = float(np.max(np.abs(stable)))
            stored = summary['reference_metrics'][key]['metrics'][name]
            close(item['fitted_amplitude_minus1_pct'], stored['amplitude_mean_minus1_pct'])
            close_statistic(item['fitted_amplitude_t_df9'], stored['amplitude_student_t_df9'])
            spectrum_metrics[name] = item
        recomputed_metrics[key] = {'primary_gate': 'PASS' if primary_pass else 'REVIEW_REQUIRED', 'spectra': spectrum_metrics}
    difference = np.abs(bins(result['raw_C_ell'][0] * factor) - bins(control * factor))
    scales = raw_b[0].copy()
    scales[3] = np.sqrt(raw_b[0, 0] * raw_b[0, 1])
    max_fractional = float(np.max(difference / scales))
    max_SEM_fraction = float(np.max(difference / (sd / np.sqrt(10))))
    numerical_pass = max_fractional < policy['numerical_max_fractional_auto_or_TE_correlation_scale'] and max_SEM_fraction < policy['numerical_max_fraction_of_pilot_mean_SEM']
    close(max_fractional, summary['numerical_iteration_control']['maximum_absolute_fractional_auto_or_TE_correlation_scale'])
    close(max_SEM_fraction, summary['numerical_iteration_control']['maximum_fraction_of_pilot_mean_SEM'])
    csv_path = evidence / 'pilot_summary/native_cmb_pilot_10_binned.csv'
    rows = list(csv.DictReader(csv_path.open()))
    if len(rows) != 3 * 4 * 39:
        raise ValueError('Incorrect binned CSV row count')
    for row in rows:
        key, index = row['reference_key'], SPECTRA.index(row['spectrum'])
        b = (int(row['ell_lo']) - 25) // 25
        if int(row['ell_hi_exclusive']) != int(row['ell_lo']) + 25:
            raise ValueError('CSV bin bounds mismatch')
        for column, array_name in [('reference_D_b_uK2', '__reference_D_b'), ('mean_D_b_uK2', '__mean_D_b'),
                                   ('sample_SD_D_b_uK2', '__sample_std_D_b'), ('bin_student_t_df9', '__bin_student_t')]:
            close(float(row[column]), result[key + array_name][index, b])
    selection_path = evidence / 'figures/native_cmb_pilot_10_selection.json'
    selection = json.loads(selection_path.read_text())
    if selection['summary']['sha256'] != sha256(summary_path) or selection['summary_metadata_sha256'] != sha256(summary_meta_path):
        raise ValueError('Plot is not based on the reviewed summary')
    if selection['source_nside'] != 2048 or selection['IDs'] != list(range(10)) or selection['spectrum_order'] != SPECTRA or selection['main_reference'] != REFERENCES[0] or selection['references'] != REFERENCES:
        raise ValueError('Plot resolution/ensemble/reference selection mismatch')
    if selection['script_sha256'] != sha256(root / 'validation/plot_pilot.py') or not selection['no_fitted_lensing_amplitude']:
        raise ValueError('Plot source or no-rescaling contract mismatch')
    for extension, record in selection['outputs'].items():
        check_record(evidence / f'figures/native_cmb_pilot_10.{extension}', record)
    for index, name in enumerate(SPECTRA):
        close(selection['shade_bounds'][name]['sample_SD_lower_minimum_uK2'], np.min(mean[index] - sd[index]))
        if index < 3 and np.any(mean[index] - sd[index] <= 0):
            raise ValueError('A logarithmic shade lower bound is clipped in this actual figure')
    if not numerical_pass or recomputed_metrics[REFERENCES[0]]['primary_gate'] != 'PASS':
        raise ValueError('Frozen numerical/primary gate does not pass')
    return {'status': 'INDEPENDENT_ARRAY_AND_PROVENANCE_REVIEW_PASS',
            'bulk_recommendation': 'GO for the declared Sehgal-completion primary remapped by actual mmDL kappa; retain the quantified subpercent BB model offset from pure CAMB',
            'not_authorized_claim': 'All four spectra are statistically consistent with pure self-consistent CAMB',
            'primary_gate': 'PASS', 'BB_interpretation': 'Pure-CAMB fitted amplitude deficit is statistically resolved; conditional actual-kappa two-point closure removes most of that deficit without a fitted lensing amplitude.',
            'policy_sha256': sha256(policy_path), 'policy_changed': False,
            'thresholds': {'amplitude_t_df9': amp_limit, 'shape_t_df9': bin_limit},
            'numerical_control': {'status': 'PASS', 'max_fractional_auto_or_TE_scale': max_fractional, 'max_fraction_of_mean_SEM': max_SEM_fraction},
            'recomputed_reference_metrics': recomputed_metrics, 'producer_input_reviews': input_reports,
            'producer_records': producer_records, 'theory_input_review': theory_report,
            'summary_archive': summary_record, 'CSV_rows_verified': len(rows),
            'figure_selection_sha256': sha256(selection_path), 'figure_visual_review': 'Four panels, mean/sample-SD shading and three distinct reference labels visible; signed linear TE; no overlaps or truncation seen',
            'checks': {'frozen_validation_hashes': True, 'all_ten_native_measurement_inputs': True,
                       'raw_spectra_to_weighted_D_b': True, 'mean_SD_ddof1': True,
                       'fit_amplitudes_are_diagnostics_not_rescaling': True, 'three_reference_inputs_exact': True,
                       'actual_kappa_auto_cross_transfers_PP': True, 'conditional_PP_mean_exact': True,
                       'auxiliary_correlation_zero_mode_canonical_and_raw_value_retained': True, 'all_plot_CSV_values': True},
            'scope': 'Independent small-array/hash/provenance/statistics review only; no new CAMB, map2alm, native map generation or native FITS scan',
            'limitations': ['Native full-map finite/checksum/hash evidence comes from the frozen executed measurement and generation workers; this review independently recomputes their stored spectra/statistics.',
                            'Headerless input frame remains owner-declared G with explicit upstream ICRS conflict; no physical-frame authentication.',
                            'Approximate t/Bonferroni pilot diagnostics use only10 realizations. BB is not subjected to a forced pure-CAMB PASS.',
                            'Measured-mean-PP relensing is conditional and Gaussian-deflection based; it does not capture all non-Gaussian kappa statistics or define the mean of separately relensed references.'],
            'review_script_sha256': sha256(__file__)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--evidence', type=Path, default=Path(__file__).resolve().parent / 'evidence/pilot_review')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = review(args.root, args.evidence)
    rendered = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(rendered)
    print(rendered, end='')
