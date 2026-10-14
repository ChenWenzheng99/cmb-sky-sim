#!/usr/bin/env python3
"""Independent ten-map pilot aggregation, gates and reference diagnostics."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from spectral_statistics import bin_native, evaluate, numerical_iteration_control, SPECTRA
from validate_inputs import check_record, check_theory, sha256


def summarize(root, spectra_dir, out, matched_reference=None):
    root, spectra_dir, out = Path(root), Path(spectra_dir), Path(out)
    policy_path = root / 'validation' / 'validation_policy.json'
    policy = json.loads(policy_path.read_text())
    theory_check = check_theory(root)
    theory_path = root / 'theory' / 'outputs' / 'sehgal2010_theory.npz'
    with np.load(theory_path, allow_pickle=False) as archive:
        references = {'self_consistent_CAMB': archive['lensed_raw_cl'][:1025].T.copy(),
                      'producer_bandlimit_CAMB': archive['producer_cutoff_lensed_raw_cl'].T.copy()}
    matched_record = None
    if matched_reference:
        matched_reference = Path(matched_reference)
        metadata = json.loads(matched_reference.with_suffix('.json').read_text())
        matched_record = check_record(matched_reference, metadata['artifact'])
        if metadata['primary_theory']['sha256'] != sha256(theory_path):
            raise ValueError('Matched reference uses a different primary-theory archive')
        if metadata['columns'] != list(SPECTRA) or metadata['raw_cl_units'] != 'uK_CMB^2':
            raise ValueError('Matched-reference units/columns mismatch')
        with np.load(matched_reference, allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['ell'], np.arange(1025))
            references['conditional_actual_mean_kappa_2point'] = archive['lensed_raw_cl'].T.copy()
    raw, reports, control = [], [], None
    for rid in range(10):
        path = spectra_dir / f'native_cmb_spectra_{rid:05d}.npz'
        sidecar_path = path.with_suffix('.json')
        meta = json.loads(sidecar_path.read_text())
        check_record(path, meta['output'])
        expected = {'realization_id': rid, 'source_nside': 2048, 'analysis_nside': 2048,
                    'lmax': 1024, 'iter': 3, 'pol': True, 'spectra_order': list(SPECTRA),
                    'use_weights': False, 'use_pixel_weights': False, 'beam': False,
                    'mask': False, 'pixel_window': False, 'downgrade': False,
                    'monopole_dipole_removal': False, 'temperature_unit': 'uK_CMB'}
        for key, value in expected.items():
            if meta.get(key) != value:
                raise ValueError(f'ID{rid} measurement contract mismatch: {key}')
        if meta.get('FITS_checksums_verified') is not True:
            raise ValueError('Native FITS checksum verification missing')
        producer = meta['producer']
        if producer['realization_id'] != rid or producer['output']['sha256'] != meta['source']['sha256']:
            raise ValueError('Measured map/producer realization identity mismatch')
        if producer['theory_archive']['sha256'] != sha256(theory_path):
            raise ValueError('Pilot map uses a different theoretical primary')
        for stats in meta['statistics']:
            if stats['pixels'] != 50331648 or stats['finite_count'] != 50331648 or stats['unseen_count'] != 0:
                raise ValueError('Native full-array finite check missing/failed')
        with np.load(path, allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['ell'], np.arange(1025))
            np.testing.assert_array_equal(archive['spectra'], SPECTRA)
            spectrum = archive['C_ell'].copy()
            if rid == 0:
                control = archive['control1024_iter1_C_ell'].copy()
            np.testing.assert_array_equal(archive['C_ell_six'][:4], spectrum)
            factors = np.arange(1025) * (np.arange(1025) + 1.) / (2 * np.pi)
            np.testing.assert_array_equal(archive['D_ell'], spectrum * factors)
            np.testing.assert_array_equal(archive['bin_lo'], np.arange(25, 1000, 25))
            np.testing.assert_array_equal(archive['bin_hi_exclusive'], np.arange(50, 1025, 25))
            np.testing.assert_allclose(archive['C_b'], bin_native(spectrum), rtol=1e-14, atol=0)
            np.testing.assert_allclose(archive['D_b'], bin_native(spectrum * factors), rtol=1e-14, atol=0)
            ell = np.arange(1025)
            centers = np.array([np.average(ell[lo:lo + 25], weights=2 * ell[lo:lo + 25] + 1) for lo in range(25, 1000, 25)])
            np.testing.assert_array_equal(archive['ell_bin'], centers)
        raw.append(spectrum)
        reports.append({'realization_id': rid, 'archive_sha256': sha256(path),
                        'measurement_metadata_sha256': sha256(sidecar_path),
                        'map_sha256': meta['source']['sha256'], 'seed': producer['seed'],
                        'source_kappa_sha256': producer['source_kappa']['sha256'],
                        'measurement_script_sha256': meta['script']['sha256']})
    for key in ['map_sha256', 'seed', 'source_kappa_sha256']:
        if len({item[key] for item in reports}) != 10:
            raise ValueError(f'Pilot has duplicate {key}')
    raw = np.array(raw)
    metrics, outputs = {}, {'ell': np.arange(1025), 'ids': np.arange(10), 'spectra': np.array(SPECTRA), 'raw_C_ell': raw}
    for key, reference in references.items():
        metrics[key], arrays = evaluate(raw, reference, policy)
        for name, value in arrays.items():
            outputs[key + '__' + name] = value
        outputs[key + '__reference_C_ell'] = reference
    main_arrays = {name: outputs['self_consistent_CAMB__' + name] for name in ['sample_std_D_b']}
    numerical = numerical_iteration_control(raw[0], control, main_arrays, policy)
    main = metrics['self_consistent_CAMB']
    report = {'scope': policy['scope'], 'status': 'PASS_NUMERICAL_AND_PRIMARY_GATE' if numerical['status'] == 'PASS' and main['primary_spectral_gate'] == 'PASS' else 'REVIEW_REQUIRED',
              'BB_scientific_agreement_status': 'REQUIRES_EXPLICIT_REVIEW_OF_REPORTED_BB_METRICS',
              'policy': policy, 'policy_sha256': sha256(policy_path), 'theory_input_validation': theory_check,
              'measurement_provenance': reports, 'numerical_iteration_control': numerical,
              'reference_metrics': metrics, 'matched_reference': matched_record,
              'matched_reference_available': 'conditional_actual_mean_kappa_2point' in references,
              'main_reference': 'self_consistent_CAMB',
              'reference_labels': {'self_consistent_CAMB': 'Sehgal2010 parameter completion: self-consistent CAMB lensed CMB',
                                   'producer_bandlimit_CAMB': 'Same CAMB primary, cut at5999; CAMB phi cut at6143',
                                   'conditional_actual_mean_kappa_2point': 'Same bandlimited primary; conditional actual ten-map mean kappa two-point Gaussian-deflection reference'},
              'mean_and_shading': 'Arithmetic mean of10; sample SD ddof1. SEM appears only in named statistical tests.',
              'script_sha256': sha256(__file__), 'statistics_script_sha256': sha256(Path(__file__).with_name('spectral_statistics.py'))}
    out.mkdir(parents=True, exist_ok=True)
    destination = out / 'native_cmb_pilot_10.npz'
    if destination.exists():
        raise FileExistsError(f'Use a fresh pilot-summary directory: {destination}')
    np.savez_compressed(destination, **outputs)
    report['artifact'] = {'path': str(destination.resolve()), 'bytes': destination.stat().st_size, 'sha256': sha256(destination)}
    report_path = destination.with_suffix('.json')
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    csv_path = out / 'native_cmb_pilot_10_binned.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['reference_key', 'spectrum', 'ell_lo', 'ell_hi_exclusive', 'reference_D_b_uK2',
                         'mean_D_b_uK2', 'sample_SD_D_b_uK2', 'bin_student_t_df9'])
        for reference_key in references:
            for index, label in enumerate(SPECTRA):
                for b, lower in enumerate(range(25, 1000, 25)):
                    writer.writerow([reference_key, label, lower, lower + 25,
                                     outputs[reference_key + '__reference_D_b'][index, b],
                                     outputs[reference_key + '__mean_D_b'][index, b],
                                     outputs[reference_key + '__sample_std_D_b'][index, b],
                                     outputs[reference_key + '__bin_student_t'][index, b]])
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--spectra-dir', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--matched-reference', type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(args.root, args.spectra_dir, args.out, args.matched_reference), indent=2, allow_nan=False))
