#!/usr/bin/env bash
# Source this file; the extension always uses its own package/runtime prefix.
CMB_RELEASE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ -L "${CMB_RELEASE_ROOT}/runtime" ]]; then
  printf 'Refusing a symbolic-link runtime; extract a dedicated CMB runtime.\n' >&2
  return 1
fi
CMB_RUNTIME="${CMB_RELEASE_ROOT}/runtime"
if [[ ! -x "${CMB_RUNTIME}/bin/python" ]]; then
  printf 'Missing CMB runtime. Run %s/setup_cmb_env.sh first.\n' "${CMB_RELEASE_ROOT}" >&2
  return 1
fi
CMB_RUNTIME=$(cd -- "${CMB_RUNTIME}" && pwd -P)
for CMB_PREFIX_MARKER in .ali_sky_unpacked_prefix .cmb_runtime_ready_prefix; do
  if [[ -f "${CMB_RUNTIME}/${CMB_PREFIX_MARKER}" ]] &&
     [[ "$(cat -- "${CMB_RUNTIME}/${CMB_PREFIX_MARKER}")" != "${CMB_RUNTIME}" ]]; then
    printf 'The unpacked CMB runtime was moved. Install a fresh archive at the final prefix.\n' >&2
    return 1
  fi
done
export CMB_RELEASE_ROOT CMB_RUNTIME
export PATH="${CMB_RUNTIME}/bin:${PATH}"
export PYTHONPATH="${CMB_RELEASE_ROOT}"
export PYTHONNOUSERSITE=1
unset PYTHONHOME
export LD_LIBRARY_PATH="${CMB_RUNTIME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export XDG_CACHE_HOME="${CMB_RELEASE_ROOT}/.cache"
export MPLCONFIGDIR="${XDG_CACHE_HOME}/matplotlib"
export NUMBA_CACHE_DIR="${XDG_CACHE_HOME}/numba"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p -- "${XDG_CACHE_HOME}" "${MPLCONFIGDIR}" "${NUMBA_CACHE_DIR}"
