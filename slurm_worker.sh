#!/usr/bin/env bash
# Slurm copies this script to its spool; the source directory is an argument.
set -euo pipefail
if [[ $# -lt 3 ]]; then
  printf 'Usage: slurm_worker.sh PACKAGE_ROOT noise|observe --config CONFIG [CLI arguments]\n' >&2
  exit 2
fi
ALI_SKY_PACKAGE_ROOT=$1
ALI_SKY_JOB_MODE=$2
shift 2
case "${ALI_SKY_JOB_MODE}" in
  noise|observe) ;;
  *) printf 'Expected noise or observe worker mode.\n' >&2; exit 2 ;;
esac
if [[ ! -x "${ALI_SKY_PACKAGE_ROOT}/run_${ALI_SKY_JOB_MODE}.sh" ]]; then
  printf 'Source launcher was not found in: %s\n' "${ALI_SKY_PACKAGE_ROOT}" >&2
  exit 1
fi
cd -- "${ALI_SKY_PACKAGE_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export ALI_SKY_MPI_RANKS=${ALI_SKY_MPI_RANKS:-${SLURM_NTASKS:-5}}
printf 'Job %s, host %s, allocated tasks %s, coordinator %s\n' \
  "${SLURM_JOB_ID:-interactive}" "$(hostname)" "${SLURM_NTASKS:-unknown}" "${ALI_SKY_JOB_MODE}"
exec "${ALI_SKY_PACKAGE_ROOT}/run_${ALI_SKY_JOB_MODE}.sh" "$@"
