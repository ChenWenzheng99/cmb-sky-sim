#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --job-name=sehgal2010_inputs
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/validate_theory_%j.log
set -euo pipefail
CMB_RELEASE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010
export PATH="${CMB_RELEASE_ROOT}/runtime/bin:${PATH}"
export LD_LIBRARY_PATH="${CMB_RELEASE_ROOT}/runtime/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${CMB_RELEASE_ROOT}" PYTHONNOUSERSITE=1
unset PYTHONHOME
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${CMB_RELEASE_ROOT}/.cache/matplotlib"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/validation/validate_inputs.py" --root "${CMB_RELEASE_ROOT}" --report "${CMB_RELEASE_ROOT}/evidence/theory_independent_validation.json"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -m unittest discover -s "${CMB_RELEASE_ROOT}/production" -p 'test_*.py' -v
"${CMB_RELEASE_ROOT}/runtime/bin/python" -m unittest discover -s "${CMB_RELEASE_ROOT}/validation" -p 'test_*.py' -v
