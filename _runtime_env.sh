#!/usr/bin/env bash
# Shared environment for the source package and its unpacked runtime.
ALI_SKY_PACKAGE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ALI_SKY_RUNTIME=${ALI_SKY_RUNTIME:-"${ALI_SKY_PACKAGE_ROOT}/runtime"}
if [[ ! -x "${ALI_SKY_RUNTIME}/bin/python" ]]; then
  printf 'Missing runtime: %s\nRun setup_env.sh first.\n' "${ALI_SKY_RUNTIME}" >&2
  return 1
fi
export ALI_SKY_PACKAGE_ROOT ALI_SKY_RUNTIME
export PATH="${ALI_SKY_RUNTIME}/bin:${PATH}"
export PYTHONPATH="${ALI_SKY_PACKAGE_ROOT}"
export PYTHONNOUSERSITE=1
unset PYTHONHOME
export LD_LIBRARY_PATH="${ALI_SKY_RUNTIME}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYSM_LOCAL_DATA=${PYSM_LOCAL_DATA:-"${ALI_SKY_PACKAGE_ROOT}/template_cache"}
export XDG_CACHE_HOME=${XDG_CACHE_HOME:-"${ALI_SKY_PACKAGE_ROOT}/.cache"}
export NUMBA_CACHE_DIR=${NUMBA_CACHE_DIR:-"${XDG_CACHE_HOME}/numba"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-"${XDG_CACHE_HOME}/matplotlib"}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMBA_NUM_THREADS=${NUMBA_NUM_THREADS:-1}
export ALI_SKY_MPIEXEC="${ALI_SKY_RUNTIME}/bin/mpiexec"
export ALI_SKY_MPI_RANKS=${ALI_SKY_MPI_RANKS:-5}
mkdir -p -- "${XDG_CACHE_HOME}" "${NUMBA_CACHE_DIR}" "${MPLCONFIGDIR}"
