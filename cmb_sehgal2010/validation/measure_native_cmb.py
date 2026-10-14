#!/usr/bin/env python3
"""Native NSIDE2048 T/Q/U -> TT/EE/BB/TE spectrum, no downgrade or transfer."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time

from astropy.io import fits
import healpy as hp
import numpy as np


NSIDE, LMAX, ITER = 2048, 1024, 3


def file_record(path):
    path = Path(path).resolve()
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'File changed while hashing: {path}')
    return {'path': str(path), 'bytes': before.st_size, 'mtime_ns': before.st_mtime_ns,
            'sha256': h.hexdigest()}


def native_stats(value):
    result = []
    for label, field in zip(('T', 'Q', 'U'), value):
        finite = unseen = 0
        minimum, maximum = np.inf, -np.inf
        for start in range(0, field.size, 1048576):
            block = field[start:start + 1048576]
            finite += int(np.isfinite(block).sum())
            unseen += int(hp.mask_bad(block).sum())
            minimum, maximum = min(minimum, float(block.min())), max(maximum, float(block.max()))
        if finite != field.size or unseen:
            raise ValueError(f'{label} has nonfinite or UNSEEN pixels')
        result.append({'field': label, 'pixels': field.size, 'finite_count': finite,
                       'unseen_count': unseen, 'minimum': minimum, 'maximum': maximum})
    return result


def harmonic_spectra(value, iteration):
    alm = hp.map2alm(value, lmax=LMAX, iter=iteration, pol=True,
                     use_weights=False, use_pixel_weights=False)
    six = np.asarray(hp.alm2cl(alm, lmax=LMAX))
    # healpy's diagonal ordering is TT,EE,BB,TE,EB,TB.
    if six.shape != (6, LMAX + 1) or not np.isfinite(six).all():
        raise ValueError('Unexpected native harmonic spectrum')
    del alm
    return six


def main(args):
    start = time.monotonic()
    source = file_record(args.map)
    header = fits.getheader(args.map, 1)
    if header.get('NSIDE') != NSIDE or header.get('ORDERING') != 'RING' or header.get('TFIELDS') != 3:
        raise ValueError('Expected native NSIDE2048 RING three-field FITS')
    if header.get('COORDSYS') not in ('G', 'C') or header.get('POLCCONV') != 'COSMO':
        raise ValueError('New release output must state G/C and COSMO')
    if any(header.get(f'TUNIT{i}') != 'uK_CMB' for i in (1, 2, 3)):
        raise ValueError('Expected explicit thermodynamic uK_CMB')
    if [header.get(f'TTYPE{i}') for i in (1, 2, 3)] != ['T', 'Q', 'U'] or any(not header.get(f'TFORM{i}', '').endswith('D') for i in (1, 2, 3)):
        raise ValueError('Expected float64 T/Q/U columns')
    with fits.open(args.map, memmap=True) as hdus:
        checksum_valid = all('CHECKSUM' in hdu.header and hdu.verify_checksum() == 1 and hdu.verify_datasum() == 1 for hdu in hdus)
    if not checksum_valid:
        raise ValueError('New release FITS checksums are missing or invalid')
    producer_path = args.producer_metadata or args.map.with_suffix('.json')
    producer = json.loads(producer_path.read_text())
    if args.id not in range(10) or producer.get('realization_id') != args.id:
        raise ValueError('Pilot realization ID must match the producer sidecar and lie in 0..9')
    if producer['coordinate_provenance']['output_coordinate'] != header['COORDSYS']:
        raise ValueError('Producer/FITS coordinate declaration mismatch')
    # The producer sidecar must tie the measured output to its actual bytes.
    if producer.get('output', {}).get('sha256') != source['sha256'] or producer['output'].get('bytes') != source['bytes']:
        raise ValueError('Producer/output identity mismatch')
    value = hp.read_map(str(args.map), field=(0, 1, 2), dtype=np.float64, nest=False)
    if value.shape != (3, 50331648):
        raise ValueError('Native map shape mismatch')
    statistics = native_stats(value)
    six = harmonic_spectra(value, ITER)
    ell = np.arange(LMAX + 1)
    factor = ell * (ell + 1.) / (2 * np.pi)
    lo = np.arange(25, 1000, 25)
    hi = lo + 25
    weights = 2 * ell + 1
    def bins(arr):
        return np.array([np.sum(arr[:, lower:upper] * weights[lower:upper], axis=1) /
                         np.sum(weights[lower:upper]) for lower, upper in zip(lo, hi)]).T
    arrays = {'ell': ell, 'C_ell': six[:4], 'D_ell': six[:4] * factor,
              'C_ell_six': six, 'bin_lo': lo, 'bin_hi_exclusive': hi,
              'ell_bin': np.array([np.average(ell[lower:upper], weights=weights[lower:upper])
                                   for lower, upper in zip(lo, hi)]), 'C_b': bins(six[:4]),
              'D_b': bins(six[:4] * factor), 'spectra': np.array(['TT', 'EE', 'BB', 'TE'])}
    if args.id == 0:
        control = harmonic_spectra(value, 1)
        arrays['control1024_iter1_C_ell'] = control[:4]
        arrays['control1024_iter1_C_ell_six'] = control
    after = Path(args.map).stat()
    if (source['bytes'], source['mtime_ns']) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Source changed during native transform')
    args.out.mkdir(parents=True, exist_ok=True)
    destination = args.out / f'native_cmb_spectra_{args.id:05d}.npz'
    if destination.exists():
        raise FileExistsError(f'Use a fresh validation output directory: {destination}')
    temporary = destination.with_name(destination.name + f'.tmp.{os.getpid()}.npz')
    np.savez_compressed(temporary, **arrays)
    os.link(temporary, destination)
    temporary.unlink()
    report = {'realization_id': args.id, 'source': source, 'statistics': statistics,
              'producer_metadata': file_record(producer_path), 'producer': producer,
              'source_nside': NSIDE, 'analysis_nside': NSIDE, 'lmax': LMAX, 'iter': ITER,
              'pol': True, 'alm_order': ['T', 'E', 'B'], 'spectra_order': ['TT', 'EE', 'BB', 'TE'],
              'use_weights': False, 'use_pixel_weights': False,
              'beam': False, 'mask': False, 'pixel_window': False, 'downgrade': False,
              'monopole_dipole_removal': False, 'temperature_unit': 'uK_CMB',
              'input_output_coordinate': header['COORDSYS'],
              'FITS_checksums_verified': checksum_valid,
              'coordinate_certification': 'Recorded producer/FITS declaration; full-sky spectra do not authenticate a headerless source frame',
              'display': {'bin_quantity': 'D_ell=ell*(ell+1)*C_ell/(2*pi)',
                          'weight': '2*ell+1', 'bin_width': 25, 'ell_interval': '[25,1000)'},
              'ID0_control': 'same native map, lmax1024, iter1' if args.id == 0 else None,
              'script': file_record(Path(__file__)), 'output': file_record(destination),
              'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'healpy': hp.__version__},
              'slurm_job_id': os.environ.get('SLURM_JOB_ID'), 'host': platform.node(),
              'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
              'elapsed_seconds': time.monotonic() - start}
    destination.with_suffix('.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': 'MEASURED', 'id': args.id, 'output': str(destination)}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--id', required=True, type=int)
    p.add_argument('--map', required=True, type=Path)
    p.add_argument('--producer-metadata', type=Path)
    p.add_argument('--out', required=True, type=Path)
    main(p.parse_args())
