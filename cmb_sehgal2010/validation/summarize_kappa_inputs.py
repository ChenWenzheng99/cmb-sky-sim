#!/usr/bin/env python3
"""Audit ten existing producer input spectra and expose their mean phi PP.

No source map is loaded and no spherical harmonic transform is performed.
The small output is the explicit input for a conditional continuum reference.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from validate_inputs import check_producer_input_spectra, sha256


def summarize(producer_dir, theory, out):
    producer_dir, theory, out = Path(producer_dir), Path(theory), Path(out)
    reports, inputs, raw_pp, seeds, source_hashes = [], [], [], [], []
    for rid in range(10):
        stem = f'lensedcmb_new_{rid:05d}'
        candidates = [producer_dir / (stem + '.json'), producer_dir / f'{rid:05d}' / (stem + '.json')]
        existing = [path for path in candidates if path.exists()]
        if len(existing) != 1:
            raise ValueError(f'Expected exactly one producer sidecar for ID{rid}: {candidates}')
        metadata_path = existing[0]
        meta = json.loads(metadata_path.read_text())
        if meta['realization_id'] != rid:
            raise ValueError('Producer realization ID mismatch')
        report, arrays, pp = check_producer_input_spectra(metadata_path.with_name(stem + '.input_spectra.npz'), meta, theory)
        reports.append(report)
        inputs.append(arrays)
        raw_pp.append(pp)
        seeds.append(meta['seed'])
        source_hashes.append(meta['source_kappa']['sha256'])
    if len(set(seeds)) != 10 or len(set(source_hashes)) != 10:
        raise ValueError('Pilot seeds or kappa input byte identities are not distinct')
    raw_pp = np.array(raw_pp)
    arrays = {'ell': np.arange(6144, dtype=np.int64), 'ids': np.arange(10, dtype=np.int64),
              'raw_phi_PP': raw_pp, 'mean_phi_PP': raw_pp.mean(0),
              'sample_std_phi_PP': raw_pp.std(0, ddof=1),
              'raw_kappa_C_ell': np.array([item['kappa_C_ell'] for item in inputs]),
              'raw_phi_C_ell': np.array([item['phi_C_ell'] for item in inputs]),
              'raw_deflection_C_ell': np.array([item['deflection_C_ell'] for item in inputs])}
    out.mkdir(parents=True, exist_ok=True)
    destination = out / 'producer_kappa_ensemble_10.npz'
    if destination.exists():
        raise FileExistsError(f'Choose a fresh output directory: {destination}')
    np.savez_compressed(destination, **arrays)
    report = {'status': 'PASS', 'ids': list(range(10)), 'input_reports': reports,
              'seeds': seeds, 'source_kappa_sha256': source_hashes,
              'theory': {'path': str(theory.resolve()), 'sha256': sha256(theory)},
              'artifact': {'path': str(destination.resolve()), 'sha256': sha256(destination), 'bytes': destination.stat().st_size},
              'script_sha256': sha256(__file__), 'sample_standard_deviation_ddof': 1,
              'PP_definition': '[L(L+1)]^2*C_phi/(2*pi); dimensionless, zero-based ell',
              'scope': 'Existing native2048 producer-alm spectra only; no further kappa transform',
              'limitations': ['Actual ten-map PP is conditional input, not an independent prediction.',
                              'Mean PP relensing differs from mean separately relensed predictions.',
                              'No statement that a two-point Gaussian-deflection reference reproduces connected non-Gaussian kappa statistics.']}
    destination.with_suffix('.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-dir', required=True, type=Path)
    parser.add_argument('--theory', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(args.producer_dir, args.theory, args.out), indent=2, allow_nan=False))
