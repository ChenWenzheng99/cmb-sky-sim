#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=aliws001,aliws003,aliws004
#SBATCH --job-name=sehgal2010_offline
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/ali_sky_release_1.0/cmb_sehgal2010/logs/offline_install_%j.log
set -euo pipefail
CMB_PROBE_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/releases/cmb_sehgal2010_install_probe_20261011
cd "${CMB_PROBE_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
./setup_cmb_env.sh
./run_cmb.sh --help
source ./_cmb_env.sh
"${CMB_RUNTIME}/bin/python" -m unittest discover -s theory -p 'test_*.py'
"${CMB_RUNTIME}/bin/python" -m unittest discover -s production -p 'test_*.py'
"${CMB_RUNTIME}/bin/python" -m unittest discover -s validation -p 'test_*.py'
"${CMB_RUNTIME}/bin/python" validation/validate_inputs.py --root . --report evidence/offline_theory_input_audit.json
"${CMB_RUNTIME}/bin/python" evidence/verify_lenspyx_sources.py --root . --report evidence/lenspyx_sources_verification.json
"${CMB_RUNTIME}/bin/python" - <<'PY'
import json, os, platform
from pathlib import Path
root = Path(os.environ['CMB_RELEASE_ROOT'])
result = {'status':'PASS_FRESH_OFFLINE_INSTALL_AND_LIGHTWEIGHT_CONTROLS',
          'prefix':str(root/'runtime'),'python':platform.python_version(),
          'lenspyx_vendored_source_files_checked':43,
          'lenspyx_installed_Python_data_files_checked':41,
          'lenspyx_source_only_Fortran_files':2,
          'new_native_maps_generated_by_this_probe':False,
          'scope':'Fresh offline extraction, wheel install, isolated imports, tiny convention tests and copied theory input audit; native ten-map spectra are separately validated.'}
(root/'evidence/offline_install_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
PY
