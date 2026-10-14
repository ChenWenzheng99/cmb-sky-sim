#!/usr/bin/env bash
set -euo pipefail
CMB_RELEASE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ $# -gt 1 ]]; then
  printf 'Usage: %s [runtime_base.tar.gz]\n' "$0" >&2
  exit 2
fi
if [[ ${1:-} == --help || ${1:-} == -h ]]; then
  printf 'Usage: %s [runtime_base.tar.gz]\nInstalls offline into this package/runtime at its final path.\n' "$0"
  exit 0
fi
CMB_BASE_ARCHIVE=${1:-"${CMB_RELEASE_ROOT}/runtime_base.tar.gz"}
CMB_RUNTIME="${CMB_RELEASE_ROOT}/runtime"
CMB_BASE_SHA=fa263b88f8f101e71d0f7de31cc68593202fa5d52ee7d5dc55b279d7dd94c577
if [[ ! -f "${CMB_BASE_ARCHIVE}" ]]; then
  printf 'Missing immutable base archive: %s\n' "${CMB_BASE_ARCHIVE}" >&2
  exit 1
fi
read -r CMB_ACTUAL_SHA _ < <(sha256sum -- "${CMB_BASE_ARCHIVE}")
if [[ ${CMB_ACTUAL_SHA} != "${CMB_BASE_SHA}" ]]; then
  printf 'Base runtime SHA256 mismatch; no installation performed.\n' >&2
  exit 1
fi
if [[ -L "${CMB_RUNTIME}" ]]; then
  printf 'Refusing a symbolic-link runtime; the original release must remain separate.\n' >&2
  exit 1
fi
if command -v flock >/dev/null 2>&1; then
  exec 9>"${CMB_RELEASE_ROOT}/.cmb_setup.lock"
  if ! flock -n 9; then
    printf 'Another setup process holds this package lock.\n' >&2
    exit 1
  fi
fi
if [[ ! -x "${CMB_RUNTIME}/bin/python" ]]; then
  if [[ -d "${CMB_RUNTIME}" ]] && [[ -n "$(ls -A -- "${CMB_RUNTIME}")" ]]; then
    printf 'Runtime destination is not empty: %s\n' "${CMB_RUNTIME}" >&2
    exit 1
  fi
  mkdir -p -- "${CMB_RUNTIME}"
  tar -xzf "${CMB_BASE_ARCHIVE}" -C "${CMB_RUNTIME}"
elif [[ ! -f "${CMB_RUNTIME}/.ali_sky_unpacked_prefix" ]]; then
  printf 'Existing runtime has no verified unpacked-prefix marker. Use a fresh package directory.\n' >&2
  exit 1
fi
source "${CMB_RELEASE_ROOT}/_cmb_env.sh"
CMB_PREFIX_MARKER="${CMB_RUNTIME}/.ali_sky_unpacked_prefix"
if [[ ! -f "${CMB_PREFIX_MARKER}" ]]; then
  "${CMB_RUNTIME}/bin/python" "${CMB_RUNTIME}/bin/conda-unpack"
  printf '%s\n' "${CMB_RUNTIME}" >"${CMB_PREFIX_MARKER}"
fi
"${CMB_RUNTIME}/bin/python" - <<'PY'
import hashlib, json, os
from pathlib import Path
root = Path(os.environ['CMB_RELEASE_ROOT'])
manifest = root / 'evidence/offline_wheels_manifest.json'
records = json.loads(manifest.read_text())['files']
expected = sorted(p.name for p in (root / 'wheels').glob('*.whl'))
assert len(records) == 7, 'Expected seven audited offline wheels'
assert sorted(record['path'] for record in records) == expected, 'Wheel inventory differs from manifest'
for record in records:
    path = root / 'wheels' / record['path']
    assert path.is_file() and path.stat().st_size == record['bytes'], f'Wheel size mismatch: {path.name}'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256'], f'Wheel SHA256 mismatch: {path.name}'
print('Seven offline wheel hashes verified')
PY
rm -f -- "${CMB_RUNTIME}/.cmb_runtime_ready_prefix"
"${CMB_RUNTIME}/bin/python" -m pip install --no-index --no-deps --force-reinstall \
  --find-links "${CMB_RELEASE_ROOT}/wheels" -r "${CMB_RELEASE_ROOT}/requirements-extra.txt"
"${CMB_RUNTIME}/bin/python" -m pip check
"${CMB_RUNTIME}/bin/python" - <<'PY'
import importlib.metadata as metadata
import json, platform
import camb, lenspyx, healpy, numpy, scipy, ducc0, contourpy, sympy, mpmath
from mpi4py import MPI
expected = {'numpy':'1.26.4', 'camb':'1.6.1', 'ducc0':'0.38.0',
            'contourpy':'1.3.3', 'sympy':'1.14.0', 'mpmath':'1.3.0', 'lenspyx':'2.0.0'}
versions = {name:metadata.version(name) for name in expected}
assert versions == expected, f'Installed overrides differ: {versions}'
assert platform.python_version_tuple()[:2] == ('3', '12'), 'Expected Python3.12 base runtime'
print(json.dumps({'python':platform.python_version(), 'overrides':versions,
                  'healpy':metadata.version('healpy'), 'scipy':metadata.version('scipy'),
                  'MPI':MPI.get_vendor()}, indent=2))
PY
printf '%s\n' "${CMB_RUNTIME}" >"${CMB_RUNTIME}/.cmb_runtime_ready_prefix"
printf 'CMB runtime ready: %s\nKeep this unpacked runtime at its final prefix.\n' "${CMB_RUNTIME}"
