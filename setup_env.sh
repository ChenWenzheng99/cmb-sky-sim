#!/usr/bin/env bash
set -euo pipefail
ALI_SKY_PACKAGE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ALI_SKY_RUNTIME=${ALI_SKY_RUNTIME:-"${ALI_SKY_PACKAGE_ROOT}/runtime"}
ALI_SKY_RUNTIME_ARCHIVE=${1:-"${ALI_SKY_PACKAGE_ROOT}/runtime.tar.gz"}
if [[ $# -gt 1 ]]; then
  printf 'Usage: %s [runtime.tar.gz]\n' "$0" >&2
  exit 2
fi
if [[ ! -x "${ALI_SKY_RUNTIME}/bin/python" ]]; then
  if [[ ! -f "${ALI_SKY_RUNTIME_ARCHIVE}" ]]; then
    printf 'Runtime archive not found: %s\n' "${ALI_SKY_RUNTIME_ARCHIVE}" >&2
    exit 1
  fi
  if [[ -d "${ALI_SKY_RUNTIME}" ]] && [[ -n "$(ls -A -- "${ALI_SKY_RUNTIME}")" ]]; then
    printf 'Runtime destination is not empty: %s\n' "${ALI_SKY_RUNTIME}" >&2
    exit 1
  fi
  mkdir -p -- "${ALI_SKY_RUNTIME}"
  tar -xzf "${ALI_SKY_RUNTIME_ARCHIVE}" -C "${ALI_SKY_RUNTIME}"
fi
ALI_SKY_RUNTIME=$(cd -- "${ALI_SKY_RUNTIME}" && pwd -P)
export ALI_SKY_RUNTIME
source "${ALI_SKY_PACKAGE_ROOT}/_runtime_env.sh"
ALI_SKY_PREFIX_MARKER="${ALI_SKY_RUNTIME}/.ali_sky_unpacked_prefix"
if [[ -f "${ALI_SKY_PREFIX_MARKER}" ]]; then
  if [[ "$(cat -- "${ALI_SKY_PREFIX_MARKER}")" != "${ALI_SKY_RUNTIME}" ]]; then
    printf 'This unpacked runtime was moved. Extract a fresh archive at the final path.\n' >&2
    exit 1
  fi
else
  "${ALI_SKY_RUNTIME}/bin/python" "${ALI_SKY_RUNTIME}/bin/conda-unpack"
  printf '%s\n' "${ALI_SKY_RUNTIME}" >"${ALI_SKY_PREFIX_MARKER}"
fi
"${ALI_SKY_RUNTIME}/bin/python" -m pip check
"${ALI_SKY_RUNTIME}/bin/python" - <<'PY'
import importlib.metadata as metadata
from mpi4py import MPI
import numpy, scipy, healpy, astropy, pysm3, numba, h5py
for name in ['numpy', 'scipy', 'healpy', 'astropy', 'pysm3', 'numba', 'mpi4py']:
    print(f'{name}={metadata.version(name)}')
print('MPI:', MPI.get_vendor())
PY
printf 'Runtime ready: %s\nPySM templates: %s\n' "${ALI_SKY_RUNTIME}" "${PYSM_LOCAL_DATA}"
