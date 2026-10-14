#!/usr/bin/env bash
# Submit one serial coordinator; the CLI launches matching MPI workers.
set -euo pipefail
ALI_SKY_PACKAGE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
ALI_SKY_JOB_MODE=${1:-}
shift
case "${ALI_SKY_JOB_MODE}" in
  noise) ALI_SKY_JOB_DEFAULT_MEM=64G; ALI_SKY_JOB_DEFAULT_TIME=12:00:00 ;;
  observe) ALI_SKY_JOB_DEFAULT_MEM=96G; ALI_SKY_JOB_DEFAULT_TIME=24:00:00 ;;
  *) printf 'Expected noise or observe submission mode.\n' >&2; exit 2 ;;
esac
if [[ $# -eq 0 ]] || [[ "${1:-}" == --help ]] || [[ "${1:-}" == -h ]]; then
  printf 'Usage: submit_%s.sh --config CONFIG.json [CLI arguments]\n' "${ALI_SKY_JOB_MODE}"
  cat <<'HELP'
Pass-through CLI options: --n-realizations, --start-id, --seed,
--output-dir, --mpi-ranks. Configuration paths resolve before submission.
Resource overrides: ALI_SKY_JOB_RANKS (default 5), ALI_SKY_JOB_MEM,
ALI_SKY_JOB_TIME, ALI_SKY_JOB_PARTITION (ali), ALI_SKY_JOB_QOS (regular),
ALI_SKY_JOB_LOG_DIR (package/logs), ALI_SKY_JOB_EXCLUDE (additional nodes).
aliws001,aliws003,aliws004 are always excluded. One node and one CPU per
task are requested.
Keep JSON/CLI mpi_ranks consistent with the requested task allocation.
HELP
  if [[ $# -eq 0 ]]; then exit 2; fi
  exit 0
fi
ALI_SKY_CLI_ARGS=()
ALI_SKY_CONFIG_SEEN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 ]]; then printf '%s\n' '--config requires a path.' >&2; exit 2; fi
      ALI_SKY_CONFIG_PATH=$2
      shift 2
      ;;
    --config=*)
      ALI_SKY_CONFIG_PATH=${1#--config=}
      shift
      ;;
    *) ALI_SKY_CLI_ARGS+=("$1"); shift; continue ;;
  esac
  if [[ ${ALI_SKY_CONFIG_SEEN} -ne 0 ]]; then
    printf 'Specify --config once.\n' >&2; exit 2
  fi
  if [[ ! -f "${ALI_SKY_CONFIG_PATH}" ]] || [[ ! -r "${ALI_SKY_CONFIG_PATH}" ]]; then
    printf 'Configuration is not a readable file: %s\n' "${ALI_SKY_CONFIG_PATH}" >&2; exit 2
  fi
  ALI_SKY_CONFIG_PATH=$(cd -- "$(dirname -- "${ALI_SKY_CONFIG_PATH}")" && pwd -P)/$(basename -- "${ALI_SKY_CONFIG_PATH}")
  ALI_SKY_CLI_ARGS+=(--config "${ALI_SKY_CONFIG_PATH}")
  ALI_SKY_CONFIG_SEEN=1
done
if [[ ${ALI_SKY_CONFIG_SEEN} -eq 0 ]]; then
  printf 'A --config JSON path is required.\n' >&2; exit 2
fi
ALI_SKY_JOB_RANKS=${ALI_SKY_JOB_RANKS:-5}
if [[ ! "${ALI_SKY_JOB_RANKS}" =~ ^[1-9][0-9]*$ ]]; then
  printf 'ALI_SKY_JOB_RANKS must be a positive integer.\n' >&2; exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
  printf 'sbatch was not found. Submit from the Slurm cluster, or use run_*.sh inside an existing compute allocation.\n' >&2
  exit 1
fi
ALI_SKY_JOB_LOG_DIR=${ALI_SKY_JOB_LOG_DIR:-"${ALI_SKY_PACKAGE_ROOT}/logs"}
mkdir -p -- "${ALI_SKY_JOB_LOG_DIR}"
ALI_SKY_JOB_LOG_DIR=$(cd -- "${ALI_SKY_JOB_LOG_DIR}" && pwd -P)
ALI_SKY_SBATCH_ARGS=(
  --partition="${ALI_SKY_JOB_PARTITION:-ali}"
  --qos="${ALI_SKY_JOB_QOS:-regular}"
  --job-name="ali_sky_${ALI_SKY_JOB_MODE}"
  --nodes=1
  --ntasks="${ALI_SKY_JOB_RANKS}"
  --cpus-per-task=1
  --mem="${ALI_SKY_JOB_MEM:-${ALI_SKY_JOB_DEFAULT_MEM}}"
  --time="${ALI_SKY_JOB_TIME:-${ALI_SKY_JOB_DEFAULT_TIME}}"
  --exclude="aliws001,aliws003,aliws004${ALI_SKY_JOB_EXCLUDE:+,${ALI_SKY_JOB_EXCLUDE}}"
  --chdir="${ALI_SKY_PACKAGE_ROOT}"
  --output="${ALI_SKY_JOB_LOG_DIR}/${ALI_SKY_JOB_MODE}_%j.out"
  --error="${ALI_SKY_JOB_LOG_DIR}/${ALI_SKY_JOB_MODE}_%j.err"
  --export=ALL
)
exec sbatch "${ALI_SKY_SBATCH_ARGS[@]}" "${ALI_SKY_PACKAGE_ROOT}/slurm_worker.sh" \
  "${ALI_SKY_PACKAGE_ROOT}" "${ALI_SKY_JOB_MODE}" "${ALI_SKY_CLI_ARGS[@]}"
