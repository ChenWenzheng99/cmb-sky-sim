#!/usr/bin/env bash
set -euo pipefail
ALI_SKY_PACKAGE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
source "${ALI_SKY_PACKAGE_ROOT}/_runtime_env.sh"
exec "${ALI_SKY_RUNTIME}/bin/python" -m ali_sky.cli noise "$@"
