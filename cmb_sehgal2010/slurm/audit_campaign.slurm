#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --job-name=sehgal2010_audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/audit_%j.log
set -euo pipefail
CMB_RELEASE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010
export PATH="${CMB_RELEASE_ROOT}/runtime/bin:${PATH}"
export LD_LIBRARY_PATH="${CMB_RELEASE_ROOT}/runtime/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${CMB_RELEASE_ROOT}" PYTHONNOUSERSITE=1
unset PYTHONHOME
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/evidence/campaign_status.py" status \
  --inventory "${CMB_RELEASE_ROOT}/evidence/before_source_inventory.json" \
  --producer "${CMB_RELEASE_ROOT}/production/lensed_CMB_generation.py" \
  --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" \
  --master-seed 20261011 --owner-coordinate G \
  --output "${CMB_RELEASE_ROOT}/evidence/campaign_status_final.json" --require-complete
