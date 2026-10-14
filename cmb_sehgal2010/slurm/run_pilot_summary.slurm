#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --job-name=sehgal2010_pilot_review
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/pilot_summary_%j.log
set -euo pipefail
CMB_RELEASE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010
CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
export PATH="${CMB_RELEASE_ROOT}/runtime/bin:${PATH}"
export LD_LIBRARY_PATH="${CMB_RELEASE_ROOT}/runtime/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export PYTHONPATH="${CMB_RELEASE_ROOT}" PYTHONNOUSERSITE=1
unset PYTHONHOME
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${CMB_RELEASE_ROOT}/.cache/matplotlib"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/validation/summarize_kappa_inputs.py" --producer-dir "${CMB_DATA_ROOT}" --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" --out "${CMB_RELEASE_ROOT}/validation/kappa_pilot"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/theory/relens_primary_with_pp.py" --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" --pp-npz "${CMB_RELEASE_ROOT}/validation/kappa_pilot/producer_kappa_ensemble_10.npz" --out "${CMB_RELEASE_ROOT}/validation/matched_kappa_reference"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/validation/summarize_pilot.py" --root "${CMB_RELEASE_ROOT}" --spectra-dir "${CMB_RELEASE_ROOT}/validation/pilot_spectra" --matched-reference "${CMB_RELEASE_ROOT}/validation/matched_kappa_reference/matched_mmdl_pp_reference.npz" --out "${CMB_RELEASE_ROOT}/validation/pilot_summary"
"${CMB_RELEASE_ROOT}/runtime/bin/python" -u "${CMB_RELEASE_ROOT}/validation/plot_pilot.py" --summary "${CMB_RELEASE_ROOT}/validation/pilot_summary/native_cmb_pilot_10.npz" --out "${CMB_RELEASE_ROOT}/validation/figures"
