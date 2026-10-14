"""Small runtime/MPI relocation check; independent of archive maps."""
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import platform
import sys

import astropy
import h5py
import healpy as hp
import numba
import numpy as np
import pysm3
import scipy
from mpi4py import MPI
from pysm3.utils.data import RemoteData


comm = MPI.COMM_WORLD
prefix = Path(sys.prefix).resolve()
loaded = sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
                 if '/' in line.split()[-1]})
entry = {
    'rank': comm.rank,
    'size': comm.size,
    'hostname': platform.node(),
    'executable': str(Path(sys.executable).resolve()),
    'prefix': str(prefix),
    'MPI_vendor': MPI.get_vendor(),
    'loaded_MPI_libraries': [p for p in loaded if '/libmpi' in p or '/libopen-pal' in p],
    'old_environment_paths_loaded': [p for p in loaded if '/env_py312/' in p or '/envs/stream/' in p],
    'marker_matches_prefix': (prefix / '.ali_sky_unpacked_prefix').read_text().strip() == str(prefix),
}
gathered = comm.gather(entry, root=0)
if comm.rank == 0:
    template_dir = Path(os.environ['PYSM_LOCAL_DATA']).resolve()
    template_manifest = json.loads((template_dir / 'template_manifest.json').read_text())
    template_paths = []
    for item in template_manifest['files']:
        name = item['logical_name']
        resolved = Path(RemoteData().get(name)).resolve()
        if resolved != template_dir / name:
            raise RuntimeError(f'Template did not resolve to supplied offline cache: {name}')
        template_paths.append(str(resolved))
    cl = hp.anafast(np.ones(hp.nside2npix(8)), lmax=16, iter=3, pol=False)
    passed = (all(item['marker_matches_prefix'] and not item['old_environment_paths_loaded']
                  and item['prefix'] == str(prefix) for item in gathered)
              and np.isfinite(cl).all() and len(gathered) == comm.size)
    report = {
        'status': 'PASS' if passed else 'FAIL',
        'scope': 'Relocated Linux runtime imports, exact offline-template resolution, independent NSIDE8 scalar transform and MPI gather',
        'host': platform.node(),
        'platform': platform.platform(),
        'libc': platform.libc_ver(),
        'python': platform.python_version(),
        'versions': {name: metadata.version(name) for name in
                     ['numpy', 'scipy', 'healpy', 'astropy', 'pysm3', 'numba', 'mpi4py', 'h5py', 'setuptools']},
        'MPI_library_version': MPI.Get_library_version(),
        'rank_reports': gathered,
        'sys_path': sys.path,
        'offline_template_paths': template_paths,
        'small_transform': {'nside': 8, 'lmax': 16, 'iter': 3, 'pol': False,
                            'input': 'independently constructed all-ones scalar map',
                            'finite': bool(np.isfinite(cl).all()), 'monopole_C_ell': float(cl[0])},
    }
    Path(sys.argv[1]).write_text(json.dumps(report, indent=2) + '\n')
    print(report['status'], 'MPI ranks', comm.size, 'prefix', prefix, flush=True)
    if not passed:
        comm.Abort(1)
