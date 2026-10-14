#!/usr/bin/env python3
"""Calculate a conditional continuum reference from measured mean mmDL phi PP.

No map transform, fit, amplitude correction or old-map mutation is performed.
The explicit input PP uses [L(L+1)]**2 C_L(phi,phi)/(2*pi), zero-based L.
"""
import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import platform
import tempfile
import time


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def canonicalize_correlation_spectra(np, returned, ell):
    """Store the special CAMB zero mode separately and define raw C0=D0=0."""
    if returned.shape != (len(ell), 4) or not np.isfinite(returned).all():
        raise ValueError('Invalid correlations return shape/finite domain')
    raw_zero_mode = returned[0].copy()
    canonical_D = returned.copy()
    canonical_D[0] = 0.0
    factor = ell * (ell + 1) / (2 * np.pi)
    raw = np.divide(canonical_D, factor[:, None], out=np.zeros_like(canonical_D),
                    where=factor[:, None] > 0)
    return raw, canonical_D, raw_zero_mode


def calculate(args):
    import camb
    from camb import correlations
    import numpy as np

    out = Path(args.out).resolve()
    names = ['matched_mmdl_pp_reference.npz', 'matched_mmdl_pp_reference.json']
    if any((out / name).exists() for name in names) and not args.overwrite:
        raise FileExistsError('Reference already exists; choose a fresh output directory')

    theory_path = Path(args.theory).resolve()
    manifest_path = theory_path.parent / 'theory_manifest.json'
    manifest = json.loads(manifest_path.read_text())
    config = manifest['configuration']
    if camb.__version__ != config['expected_camb_version']:
        raise ValueError('CAMB version differs from the frozen theory configuration')
    if manifest['artifact_manifest'][theory_path.name]['sha256'] != sha256(theory_path):
        raise ValueError('Primary-theory artifact hash mismatch')
    table_lmax = config['calculation']['table_lmax']
    primary_lmax = config['producer_cutoffs']['primary_lmax']
    pp_lmax = config['producer_cutoffs']['kappa_deflection_lmax']
    output_lmax = config['calculation']['cutoff_control_output_lmax']
    with np.load(theory_path, allow_pickle=False) as data:
        ell = data['ell'].copy()
        primary_raw = data['producer_unlensed_raw_cl'].copy()
    np.testing.assert_array_equal(ell, np.arange(table_lmax + 1, dtype=np.int64))
    if primary_raw.shape != (table_lmax + 1, 4) or primary_raw.dtype != np.float64 or not np.isfinite(primary_raw).all():
        raise ValueError('Invalid explicit TT,EE,BB,TE primary array')
    if np.count_nonzero(primary_raw[primary_lmax + 1:]) or np.count_nonzero(primary_raw[:, 2]):
        raise ValueError('Primary cutoff or zero-primordial-BB contract violated')
    pp_path = Path(args.pp_npz).resolve()
    with np.load(pp_path, allow_pickle=False) as data:
        if args.pp_key not in data.files:
            raise ValueError(f'Explicit PP key not found: {args.pp_key}')
        measured_pp = data[args.pp_key].copy()
        input_ell = data[args.ell_key].copy()
    if measured_pp.shape != (pp_lmax + 1,) or measured_pp.dtype != np.float64:
        raise ValueError('Measured PP must be one float64 zero-based vector through the producer deflection cutoff')
    np.testing.assert_array_equal(input_ell, np.arange(pp_lmax + 1))
    if not np.isfinite(measured_pp).all() or np.any(measured_pp < 0) or measured_pp[0] != 0:
        raise ValueError('Measured phi PP must be finite, nonnegative, and zero at L=0')
    pp = np.zeros(table_lmax + 1, dtype=np.float64)
    pp[:pp_lmax + 1] = measured_pp
    factor = ell * (ell + 1) / (2 * np.pi)
    primary_D = primary_raw * factor[:, None]
    sampling = config['calculation']['cutoff_control_sampling_factor']
    started = time.monotonic()
    print(f'Conditional mean-PP reference, full theta, sampling factor {sampling}', flush=True)
    lensed_D = correlations.lensed_cls(primary_D, pp, lmax=table_lmax,
                                      lmax_lensed=output_lmax, theta_max=None,
                                      sampling_factor=sampling, leggaus=True, cache=True)
    lensed_raw, lensed_D, raw_zero_mode = canonicalize_correlation_spectra(np, lensed_D, ell[:output_lmax + 1])
    if lensed_raw.shape != (output_lmax + 1, 4) or not np.isfinite(lensed_raw).all():
        raise ValueError('Invalid conditional-reference output')
    np.testing.assert_allclose(lensed_D, lensed_raw * factor[:output_lmax + 1, None], rtol=1e-14, atol=0)
    arrays = {'ell': np.arange(output_lmax + 1, dtype=np.int64),
              'lensed_raw_cl': lensed_raw, 'lensed_D_ell': lensed_D,
              'input_PP_ell': ell, 'input_phi_PP': pp,
              'raw_correlation_zero_mode': raw_zero_mode}
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.matched_mmdl_pp_', dir=out) as tmp_name:
        stage = Path(tmp_name)
        np.savez_compressed(stage / names[0], **arrays)
        metadata = {
            'schema': 'sehgal2010_primary_measured_mmdl_mean_PP_reference_v1',
            'status': 'COMPLETED_CONDITIONAL_CONTINUUM_REFERENCE',
            'primary_theory': {'path': str(theory_path), 'sha256': sha256(theory_path),
                               'manifest_sha256': sha256(manifest_path), 'key': 'producer_unlensed_raw_cl'},
            'measured_PP_input': {'path': str(pp_path), 'sha256': sha256(pp_path),
                                  'key': args.pp_key, 'ell_key': args.ell_key,
                                  'normalization': '[L(L+1)]^2 C_L(phi,phi)/(2pi)',
                                  'units': 'dimensionless', 'zero_based_ell_verified': True},
            'primary_lmax': primary_lmax, 'deflection_lmax': pp_lmax,
            'table_lmax': table_lmax, 'output_lmax': output_lmax,
            'angular_range': 'full [0,pi]', 'sampling_factor': sampling,
            'columns': ['TT', 'EE', 'BB', 'TE'], 'raw_cl_units': 'uK_CMB^2',
            'correlation_zero_mode_convention': {
                'canonical_raw_C0_and_D0': 0.0, 'raw_correlation_zero_mode': raw_zero_mode.tolist(),
                'source_behavior': 'CAMB1.6.1 correlations uses lfacs[0]=1; returned ell0 is not physical D_ell',
                'positive_ell_correlation_output_unchanged': True},
            'cosmology_configuration_canonical_sha256': manifest['configuration_canonical_sha256'],
            'versions': {'camb': camb.__version__, 'numpy': np.__version__, 'python': platform.python_version()},
            'installed_correlations_sha256': sha256(inspect.getsourcefile(correlations)),
            'builder_sha256': sha256(__file__), 'host': platform.node(),
            'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'elapsed_seconds': time.monotonic() - started,
            'artifact': {'path': str(out / names[0]), 'sha256': sha256(stage / names[0]),
                         'bytes': (stage / names[0]).stat().st_size},
            'scope': 'Curved-sky Gaussian-deflection continuum reference evaluated at the selected measured mean phi PP and explicitly bandlimited Sehgal2010-completion primary.',
            'limitations': ['This conditions on the selected kappa ensemble power and is not an independent cosmological prediction.',
                            'The lensing result of mean PP is not defined as the mean of individually relensed predictions.',
                            'A power-spectrum-only Gaussian-deflection reference does not include connected higher-order statistics of non-Gaussian mmDL kappa.',
                            'No fitted rescaling, map transform, or amplitude adjustment is applied.']}
        (stage / names[1]).write_text(json.dumps(metadata, indent=2, allow_nan=False) + '\n')
        for name in names:
            os.replace(stage / name, out / name)
    print(f'COMPLETE {out / names[0]}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--theory', required=True, help='Frozen Sehgal2010 theory NPZ with adjacent manifest')
    parser.add_argument('--pp-npz', required=True, help='Producer-derived measured mean phi-PP NPZ')
    parser.add_argument('--pp-key', default='mean_phi_PP', help='Explicit zero-based float64 phi PP key')
    parser.add_argument('--ell-key', default='ell', help='Explicit input zero-based ell vector key')
    parser.add_argument('--out', required=True)
    parser.add_argument('--overwrite', action='store_true')
    calculate(parser.parse_args())


if __name__ == '__main__':
    main()
