"""JSON configuration, validation and explicit archive conventions."""
from dataclasses import asdict
import json
import math
import os
from pathlib import Path

from .noise import channel_from_dict, default_channels

FORECAST = Path('/sharefs/alicpt/users/chenwz/Testarea/FORECAST')
DEFAULT_REPOSITORY = {
    'sky_dir': str(FORECAST / 'A_NEW_WORK_ILC/mmDL_maps_2048'),
    'galactic_dir': str(FORECAST / 'A_NEW_WORK_ILC/Gal_FG'),
    'noise_dir': str(FORECAST / 'A_NEW_WORK_MAIN/sim/noise_red/Ali'),
    'mask': str(FORECAST / 'A_NEW_WORK_MAIN/sim/mask_c/AliCPT_LAT_cross_PLANCK_2048_bi_C.fits'),
    'sky_coordinate': 'G', 'noise_coordinate': 'G', 'mask_coordinate': 'C',
    'frame_basis': 'Dataset-owner declaration on 2026-10-11; historical upstream ICRS provenance is retained separately',
}
PRESETS = {'dust': 'd9', 'synchrotron': 's4', 'freefree': 'f1', 'ame': 'a1'}


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return value


def _boolean(value, name):
    if not isinstance(value, bool):
        raise ValueError(f'{name} must be true or false')
    return value


def _path(value, base):
    p = Path(value).expanduser()
    return str((base / p).resolve() if not p.is_absolute() else p.resolve())


def load_config(path, mode, overrides=None):
    path = Path(path).resolve()
    raw = json.loads(path.read_text())
    if not isinstance(raw, dict):
        raise ValueError('Configuration must be a JSON object')
    return normalize(raw, mode, path.parent, overrides)


def normalize(raw, mode, base=Path('.'), overrides=None):
    allowed = {'instrument', 'nside', 'start_id', 'n_realizations', 'ids', 'coordinate', 'seed',
               'frequencies_GHz', 'channels', 'scale_factor', 'output_dir', 'mpi_ranks', 'overwrite',
               'repository', 'noise_ids', 'cache_dir', 'galactic_presets', 'template_dir', 'noise_dir',
               'noise_pattern', 'noise_coordinate', 'add_noise', 'mask', 'mask_coordinate',
               'use_repository_mask', 'lmax', 'harmonic_iter'}
    if set(raw) - allowed:
        raise ValueError(f'Unknown configuration keys: {sorted(set(raw) - allowed)}')
    raw = dict(raw)
    if any((overrides or {}).get(k) is not None for k in ('n_realizations', 'start_id')):
        raw.pop('ids', None)
    raw = {**raw, **{k: v for k, v in (overrides or {}).items() if v is not None}}
    instrument = str(raw.get('instrument', 'LAT')).upper()
    if instrument not in ('SAT', 'LAT', 'CUSTOM'):
        raise ValueError('instrument must be SAT, LAT or CUSTOM')
    nside = _integer(raw.get('nside', 2048), 'nside', 1)
    if nside & (nside - 1):
        raise ValueError('nside must be a power of two')
    start = _integer(raw.get('start_id', 0), 'start_id')
    count = _integer(raw.get('n_realizations', 1 if mode == 'observe' else 50), 'n_realizations', 1)
    ids = raw.get('ids', list(range(start, start + count)))
    if not isinstance(ids, list) or not ids:
        raise ValueError('ids must be a nonempty list')
    ids = [_integer(x, 'realization ID') for x in ids]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate realization IDs')
    coordinate = str(raw.get('coordinate', 'G')).upper()
    if coordinate not in ('C', 'G', 'E'):
        raise ValueError('coordinate must be C, G or E')
    defaults = {c.frequency_GHz: asdict(c) for c in default_channels(instrument if instrument != 'CUSTOM' else 'LAT')}
    supplied = raw.get('channels')
    if supplied is None:
        frequencies = raw.get('frequencies_GHz', list(defaults))
        if not isinstance(frequencies, list) or not frequencies:
            raise ValueError('frequencies_GHz must be a nonempty list')
        supplied = []
        for frequency in frequencies:
            frequency = float(frequency)
            if frequency not in defaults or instrument == 'CUSTOM':
                raise ValueError(f'{frequency} GHz requires an explicit channels entry with beam/noise parameters')
            supplied.append(defaults[frequency])
    if not isinstance(supplied, list) or not supplied:
        raise ValueError('channels must be a nonempty list')
    channels, custom_noise = [], instrument == 'CUSTOM'
    for value in supplied:
        if not isinstance(value, dict) or 'frequency_GHz' not in value:
            raise ValueError('Each channel requires frequency_GHz')
        frequency = float(value['frequency_GHz'])
        baseline = defaults.get(frequency, {}) if instrument != 'CUSTOM' else {}
        merged = {**baseline, **value}
        if 'nlev_uK_arcmin' in value and 'n_red_P' not in value:
            merged['n_red_P'] = value['nlev_uK_arcmin']
        if not baseline:
            for field in ('nlev_uK_arcmin', 'beam_fwhm_arcmin'):
                if field not in value:
                    raise ValueError(f'Custom {frequency} GHz channel requires {field}')
            merged = {'N_red_T': 0., 'l_knee_T': 1., 'alpha_knee_T': 0.,
                      'n_red_P': value['nlev_uK_arcmin'], 'l_knee_P': 1.,
                      'alpha_knee_P': 0., 'include_red': False, **merged}
        beam_transfer_path = merged.pop('beam_transfer_path', None)
        c = channel_from_dict(merged)
        channels.append(asdict(c))
        noise_fields = ('nlev_uK_arcmin', 'N_red_T', 'l_knee_T', 'alpha_knee_T',
                        'n_red_P', 'l_knee_P', 'alpha_knee_P', 'include_red')
        if not baseline or any(merged.get(k) != baseline.get(k) for k in noise_fields):
            custom_noise = True
        if beam_transfer_path:
            channels[-1]['beam_transfer_path'] = _path(beam_transfer_path, base)
    if len({c['frequency_GHz'] for c in channels}) != len(channels):
        raise ValueError('Duplicate frequencies')
    result = {
        'mode': mode, 'nside': nside, 'instrument': instrument, 'coordinate': coordinate,
        'seed': _integer(raw.get('seed', 20261011), 'seed'),
        'scale_factor': _integer(raw.get('scale_factor', 1), 'scale_factor', 1),
        'ids': ids, 'channels': channels, 'overwrite': _boolean(raw.get('overwrite', False), 'overwrite'),
        'output_dir': _path(raw.get('output_dir', 'outputs/' + mode), base),
        'mpi_ranks': _integer(raw.get('mpi_ranks', int(os.environ.get('ALI_SKY_MPI_RANKS', 5))), 'mpi_ranks', 1),
    }
    if mode == 'noise':
        return result
    if mode != 'observe':
        raise ValueError('Unknown mode')
    repo = {**DEFAULT_REPOSITORY, **raw.get('repository', {})}
    for name in ('sky_dir', 'galactic_dir', 'noise_dir', 'mask'):
        if repo.get(name):
            repo[name] = _path(repo[name], base)
    for name in ('sky_coordinate', 'noise_coordinate', 'mask_coordinate'):
        if repo.get(name) not in ('C', 'G', 'E'):
            raise ValueError(f'repository.{name} must be C, G or E')
    noise_ids = raw.get('noise_ids', ids)
    if not isinstance(noise_ids, list) or len(noise_ids) != len(ids):
        raise ValueError('noise_ids must have one explicit ID for each sky ID')
    result.update({
        'repository': repo, 'custom_noise': custom_noise,
        'noise_ids': [_integer(x, 'noise ID') for x in noise_ids],
        'cache_dir': _path(raw.get('cache_dir', 'cache'), base),
        'galactic_presets': {**PRESETS, **raw.get('galactic_presets', {})},
        'template_dir': _path(raw['template_dir'], base) if raw.get('template_dir') else None,
        'noise_dir': _path(raw['noise_dir'], base) if raw.get('noise_dir') else None,
        'noise_pattern': raw.get('noise_pattern', '{instrument}/{frequency}/noise_{id:05d}.fits'),
        'noise_coordinate': raw.get('noise_coordinate'),
        'add_noise': _boolean(raw.get('add_noise', True), 'add_noise'),
        'mask': _path(raw['mask'], base) if raw.get('mask') else None,
        'mask_coordinate': raw.get('mask_coordinate'),
        'use_repository_mask': _boolean(raw.get('use_repository_mask', True), 'use_repository_mask'),
        'lmax': _integer(raw.get('lmax', 3 * nside - 1), 'lmax'),
        'harmonic_iter': _integer(raw.get('harmonic_iter', 1), 'harmonic_iter'),
    })
    if result['lmax'] > 3 * nside - 1:
        raise ValueError('lmax cannot exceed 3*nside-1')
    if result['scale_factor'] != 1 and result['add_noise'] and not result['noise_dir']:
        result['custom_noise'] = True
    for key in ('mask_coordinate', 'noise_coordinate'):
        if result[key] not in (None, 'C', 'G', 'E'):
            raise ValueError(f'{key} must be C, G or E')
    if not isinstance(result['galactic_presets'], dict) or set(result['galactic_presets']) != set(PRESETS):
        raise ValueError('galactic_presets must specify dust/synchrotron/freefree/ame')
    return result


def frequency_tag(frequency):
    return format(float(frequency), '.12g')
