#!/usr/bin/env bash
set -euo pipefail
ALI_SKY_PACKAGE_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
exec "${ALI_SKY_PACKAGE_ROOT}/_submit_job.sh" observe "$@"
