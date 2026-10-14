#!/usr/bin/env bash
set -euo pipefail
CMB_RELEASE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "${CMB_RELEASE_ROOT}/_cmb_env.sh"
if [[ ! -f "${CMB_RUNTIME}/.cmb_runtime_ready_prefix" ]]; then
  printf 'CMB offline installation is incomplete. Run setup_cmb_env.sh first.\n' >&2
  exit 1
fi
exec "${CMB_RUNTIME}/bin/python" -u "${CMB_RELEASE_ROOT}/production/lensed_CMB_generation.py" "$@"
