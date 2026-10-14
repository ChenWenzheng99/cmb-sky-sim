#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --job-name=sehgal2010_pilot
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=0-9%10
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/pilot_%A_%a.log
set -euo pipefail
CMB_RELEASE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010
export PATH="${CMB_RELEASE_ROOT}/runtime/bin:${PATH}"
export LD_LIBRARY_PATH="${CMB_RELEASE_ROOT}/runtime/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${CMB_RELEASE_ROOT}" PYTHONNOUSERSITE=1
unset PYTHONHOME
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${CMB_RELEASE_ROOT}/.cache/matplotlib"
CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
CMB_ID=${SLURM_ARRAY_TASK_ID}
CMB_ID_PADDED=$(printf '%05d' "${CMB_ID}")
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/production/lensed_CMB_generation.py" --mmdl-root "${CMB_DATA_ROOT}" --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" --ids "${CMB_ID}" --master-seed 20261011 --threads "${SLURM_CPUS_PER_TASK}" --owner-coordinate G --execute --manifest "${CMB_RELEASE_ROOT}/evidence/pilot_production_${CMB_ID_PADDED}.json"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/validation/measure_native_cmb.py" --id "${CMB_ID}" --map "${CMB_DATA_ROOT}/${CMB_ID_PADDED}/lensedcmb_new_${CMB_ID_PADDED}.fits" --out "${CMB_RELEASE_ROOT}/validation/pilot_spectra"
