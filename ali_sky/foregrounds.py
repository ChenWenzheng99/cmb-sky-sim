"""Read separate Galactic maps or evaluate and cache PySM at a new band."""
from contextlib import contextmanager
import fcntl
import gc
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import time

import numpy as np

from .config import frequency_tag
from .io import atomic_json, file_info, read_fits, sha256, write_tqu


@contextmanager
def cache_lock(path):
    """A shared-filesystem POSIX lock; independent jobs cannot race a cache."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _cache_key(component, preset, frequency, nside, template_dir):
    return {'component': component, 'preset': preset, 'frequency_GHz': float(frequency),
            'nside': nside, 'coordinate': 'G', 'unit': 'uK_CMB',
            'pysm3_version': importlib.metadata.version('pysm3'),
            'template_directory': str(Path(template_dir).resolve()) if template_dir else None}


def get_galactic_component(config, component, frequency):
    preset = config['galactic_presets'][component]
    nside = config['nside']
    if not isinstance(preset, str) or not preset:
        raise ValueError(f'Invalid PySM preset for {component}')
    name = f'fg_{component}_{preset}_{frequency_tag(frequency)}GHz_nside{nside}_G.fits'
    source = Path(config['repository']['galactic_dir']) / name
    if source.exists():
        maps, coord, info = read_fits(source, nside)
        return maps, coord, {'method': 'repository component', 'preset': preset, 'input': info}
    # The key includes PySM's installed version; never reuse a cache across model versions.
    template_dir = config.get('template_dir') or os.environ.get('PYSM_LOCAL_DATA')
    key = _cache_key(component, preset, frequency, nside, template_dir)
    digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]
    path = Path(config['cache_dir']) / 'galactic' / digest / name
    with cache_lock(path.with_suffix('.lock')):
        if path.exists():
            sidecar = path.with_suffix('.json')
            if not sidecar.exists():
                raise ValueError(f'Incomplete Galactic cache: {path}')
            saved = json.loads(sidecar.read_text())
            if saved.get('cache_key') != key or saved.get('output', {}).get('sha256') != sha256(path):
                raise ValueError(f'Galactic cache identity mismatch: {path}')
            for template in saved.get('templates', []):
                current = file_info(template['path'])
                if current['bytes'] != template['bytes'] or current['sha256'] != template['sha256']:
                    raise ValueError(f'PySM template changed after cache creation: {template["path"]}')
            maps, coord, info = read_fits(path, nside)
            return maps, coord, {'method': 'PySM cache hit', 'preset': preset, 'input': info,
                                 'cache_key': key, 'templates': saved.get('templates', [])}
        if template_dir:
            os.environ['PYSM_LOCAL_DATA'] = str(template_dir)
        import pysm3
        import pysm3.units as u
        from pysm3.utils.data import RemoteData
        accessed = {}
        original_get = RemoteData.get

        def tracked_get(self, logical_name):
            resolved = original_get(self, logical_name)
            accessed[str(logical_name)] = str(resolved)
            return resolved

        print(f'PySM: generating {component} ({preset}) at {frequency:g} GHz, NSIDE={nside}', flush=True)
        started = time.monotonic()
        RemoteData.get = tracked_get
        try:
            sky = pysm3.Sky(nside=nside, preset_strings=[preset])
            emission = sky.get_emission(frequency * u.GHz)
            maps = np.asarray(emission.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(frequency * u.GHz)).value,
                              dtype=np.float64)
            if maps.shape[0] == 1:
                expanded = np.zeros((3, maps.shape[1]), dtype=np.float64)
                expanded[0] = maps[0]
                maps = expanded
            del sky, emission
        finally:
            RemoteData.get = original_get
        templates = [{'logical_name': name, **file_info(value)} for name, value in sorted(accessed.items())]
        metadata = {'product': 'Galactic component', 'cache_key': key,
                    'component': component, 'preset': preset, 'frequency_GHz': float(frequency),
                    'templates': templates, 'generation_seconds': time.monotonic() - started,
                    'beam': 'none', 'bandpass': 'monochromatic',
                    'coordinate_operation': 'PySM native Galactic frame'}
        write_tqu(path, maps, 'G', metadata)
        info = file_info(path)
        gc.collect()
        return maps, 'G', {'method': 'PySM generated', 'preset': preset,
                           'input': info, 'cache_key': key, 'templates': templates}
