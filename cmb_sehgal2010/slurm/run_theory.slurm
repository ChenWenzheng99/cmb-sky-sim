#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --job-name=sehgal2010_theory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/theory_%j.log
set -euo pipefail
CMB_RELEASE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010
export PATH="${CMB_RELEASE_ROOT}/runtime/bin:${PATH}"
export LD_LIBRARY_PATH="${CMB_RELEASE_ROOT}/runtime/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${CMB_RELEASE_ROOT}" PYTHONNOUSERSITE=1
unset PYTHONHOME
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${CMB_RELEASE_ROOT}/.cache/matplotlib"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/theory/build_sehgal2010_theory.py" --config "${CMB_RELEASE_ROOT}/theory/sehgal2010_config.json" --out "${CMB_RELEASE_ROOT}/theory/outputs"
