"""A physical forward observation: sky beam, unsmoothed noise, then mask."""
import gc
import json
from pathlib import Path
import platform
import time

import healpy as hp
import numpy as np

from .config import frequency_tag
from .coordinates import rotate_map, rotate_mask
from .foregrounds import get_galactic_component
from .io import read_fits, read_noise, sha256, write_tqu
from .sed import scale_factor


PREFIX = {'tsz': 'tsz', 'cib': 'ir_pts', 'radio': 'rad_pts', 'ksz': 'ksz'}


def read_sky(config, component, rid, frequency=None):
    folder = Path(config['repository']['sky_dir']) / f'{rid:05d}'
    fields = 3 if component == 'cmb' else 1
    factor = 1.0
    if component == 'cmb':
        path = folder / f'lensedcmb_new_{rid:05d}.fits'
        source_frequency = None
    elif component == 'ksz':
        path = folder / f'ksz_{rid:05d}.fits'
        source_frequency = None
    else:
        # Use a true archived nominal frequency if present; never relabel a nearby band.
        exact = folder / f'{PREFIX[component]}_{frequency_tag(frequency)}ghz_{rid:05d}.fits'
        source_frequency = float(frequency) if exact.exists() else 148.0
        path = exact if exact.exists() else folder / f'{PREFIX[component]}_148ghz_{rid:05d}.fits'
        factor = scale_factor(component, frequency, source_frequency)
    maps, coord, info = read_fits(path, config['nside'], fields=fields,
                                assume_coord=config['repository']['sky_coordinate'],
                                assume_unit='uK_CMB')
    maps *= factor
    source_coord = coord
    maps = rotate_map(maps, coord, config['coordinate'])
    info.update({'source_frequency_GHz': source_frequency, 'target_frequency_GHz': frequency,
                 'SED_factor': factor, 'source_coordinate': source_coord,
                 'output_coordinate': config['coordinate'],
                 'frame_basis': config['repository'].get('frame_basis', 'explicit configuration'),
                 'historical_source_frame': 'C/ICRS in the retained provider/source-trace audit',
                 'source_frame_verified': info['coordinate_from_FITS_header'],
                 'coordinate_operation': 'none' if source_coord == config['coordinate'] else 'spin-aware pixel interpolation'})
    return maps, info


def noise_path(config, channel, noise_id):
    frequency = frequency_tag(channel['frequency_GHz'])
    if config.get('noise_dir'):
        values = {'instrument': config['instrument'], 'frequency': frequency,
                  'frequency_GHz': channel['frequency_GHz'], 'id': noise_id}
        relative = config['noise_pattern'].format(**values)
        path = Path(config['noise_dir']) / relative
        return path, config.get('noise_coordinate') or config['coordinate']
    path = Path(config['repository']['noise_dir']) / config['instrument'] / frequency / f'TQU_N_{noise_id}.npy'
    return path, config['repository']['noise_coordinate']


def read_mask(config):
    path = config.get('mask')
    declared_coord = config.get('mask_coordinate')
    if not path and config.get('use_repository_mask'):
        path = config['repository'].get('mask')
        declared_coord = config['repository'].get('mask_coordinate')
    if not path:
        return None, {'method': 'full sky', 'sky_fraction': 1.0}
    mask, coord, info = read_fits(path, config['nside'], fields=1,
                                assume_coord=declared_coord, require_unit=False)
    if np.any(mask < 0) or np.any(mask > 1):
        raise ValueError('Mask weights must lie within [0,1]')
    source_coord = coord
    mask = rotate_mask(mask, coord, config['coordinate'])
    info.update({'source_coordinate': source_coord, 'output_coordinate': config['coordinate'],
                 'coordinate_operation': 'none' if source_coord == config['coordinate'] else 'nearest-pixel mask rotation',
                 'sky_fraction': float(np.mean(mask)), 'weighted_sky_fraction': float(np.mean(mask ** 2))})
    return mask, info


def beam_sky(maps, channel, lmax, harmonic_iter=1):
    if channel.get('beam_transfer_path'):
        window = np.load(channel['beam_transfer_path'], allow_pickle=False)
        if window.ndim == 1:
            window = np.repeat(window[None, :], 3, axis=0)
        if window.shape[0] != 3 or window.shape[1] < lmax + 1 or not np.isfinite(window).all():
            raise ValueError('Custom beam must be a finite NPY B_ell vector or 3x(L+1) T/E/B array covering lmax')
        window = window[:, :lmax + 1]
        beam = {'type': 'custom transfer', 'path': channel['beam_transfer_path'],
                'sha256': sha256(channel['beam_transfer_path']), 'columns': 'T,E,B'}
    else:
        fwhm = channel['beam_fwhm_arcmin']
        if fwhm == 0:
            return maps, {'type': 'identity', 'fwhm_arcmin': 0}
        common = hp.gauss_beam(fwhm * np.pi / 10800., lmax=lmax)
        window = np.repeat(common[None, :], 3, axis=0)
        beam = {'type': 'Gaussian', 'fwhm_arcmin': fwhm}
    alms = hp.map2alm(maps, lmax=lmax, iter=harmonic_iter, pol=True)
    for field, transfer in zip(alms, window):
        hp.almxfl(field, transfer, inplace=True)
    value = hp.alm2map(alms, hp.npix2nside(maps.shape[1]), lmax=lmax, pol=True)
    del alms
    return value, {**beam, 'lmax': lmax, 'map2alm_iter': harmonic_iter,
                   'pixel_window': False, 'applied_to': 'sky only'}


def observe_one(config, channel, rid, noise_id):
    frequency = channel['frequency_GHz']
    output = Path(config['output_dir']) / frequency_tag(frequency) / f'observation_{rid:05d}_noise{noise_id:05d}_{config["coordinate"]}.fits'
    if output.exists() and not config['overwrite']:
        raise FileExistsError(f'Observation output exists: {output}; use a fresh output directory or overwrite explicitly')
    started = time.monotonic()
    print(f'Observation sky ID {rid}, noise ID {noise_id}, {frequency:g} GHz, NSIDE={config["nside"]}', flush=True)
    sky, cmb_info = read_sky(config, 'cmb', rid)
    inputs = {'cmb': cmb_info, 'galactic': {}, 'extragalactic': {}}
    for component in config['galactic_presets']:
        maps, coord, info = get_galactic_component(config, component, frequency)
        maps = rotate_map(maps, coord, config['coordinate'])
        sky += maps
        info['source_coordinate'] = coord
        info['output_coordinate'] = config['coordinate']
        info['coordinate_operation'] = 'none' if coord == config['coordinate'] else 'spin-aware pixel interpolation'
        inputs['galactic'][component] = info
        del maps
        gc.collect()
    for component in PREFIX:
        maps, info = read_sky(config, component, rid, frequency)
        sky[0] += maps
        inputs['extragalactic'][component] = info
        del maps
    sky, beam_info = beam_sky(sky, channel, config['lmax'], config['harmonic_iter'])
    if config['add_noise']:
        path, assumed = noise_path(config, channel, noise_id)
        noise, coord, info = read_noise(path, config['nside'], assume_coord=assumed)
        sky += rotate_map(noise, coord, config['coordinate'])
        info.update({'source_coordinate': coord, 'output_coordinate': config['coordinate'], 'beam_applied': False})
        inputs['noise'] = info
        del noise
    else:
        inputs['noise'] = {'method': 'disabled'}
    mask, mask_info = read_mask(config)
    if mask is not None:
        sky *= mask
        del mask
    metadata = {'product': 'forward observation T/Q/U', 'frequency_GHz': frequency,
                'sky_id': rid, 'noise_id': noise_id, 'inputs': inputs, 'mask': mask_info,
                'beam': beam_info, 'operation': 'mask * (beam * (CMB + Galactic + extragalactic) + noise)',
                'extragalactic_polarization': 'zero', 'bandpass': 'monochromatic',
                'no_downgrade': True, 'monopole_dipole_removal': False,
                'configuration': config, 'host': platform.node(),
                'versions': {'numpy': np.__version__, 'healpy': hp.__version__},
                'elapsed_seconds': time.monotonic() - started}
    return write_tqu(output, sky, config['coordinate'], metadata, config['overwrite'])


def preflight(config):
    """Resolve every expensive-task input before generating any new outputs."""
    for rid in config['ids']:
        folder = Path(config['repository']['sky_dir']) / f'{rid:05d}'
        required = [folder / f'lensedcmb_new_{rid:05d}.fits', folder / f'ksz_{rid:05d}.fits']
        for component in ('tsz', 'cib', 'radio'):
            for channel in config['channels']:
                exact = folder / f'{PREFIX[component]}_{frequency_tag(channel["frequency_GHz"])}ghz_{rid:05d}.fits'
                required.append(exact if exact.exists() else folder / f'{PREFIX[component]}_148ghz_{rid:05d}.fits')
        for path in set(required):
            if not path.exists():
                raise FileNotFoundError(f'Missing matched sky component: {path}')
    if config['add_noise'] and (config.get('noise_dir') or not config['custom_noise']):
        for channel in config['channels']:
            for noise_id in config['noise_ids']:
                path, _ = noise_path(config, channel, noise_id)
                if not path.exists():
                    raise FileNotFoundError(f'Missing noise map: {path}; supply noise_dir or custom noise parameters; IDs are never wrapped modulo 50')
    read_mask(config)  # Also check resolution, frame and weight domain before work starts.
    return {'sky_ids': config['ids'], 'noise_ids': config['noise_ids'],
            'frequencies_GHz': [c['frequency_GHz'] for c in config['channels']]}
