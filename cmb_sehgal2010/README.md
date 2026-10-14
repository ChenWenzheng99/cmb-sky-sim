# Native CMB generation with Sehgal2010 primary spectra and mmDL convergence

This extension draws unlensed CMB realizations from CAMB scalar spectra based
on the Sehgal2010 parameters, then uses lenspyx and the matched archived mmDL
convergence map to produce HEALPix NSIDE=2048 T/Q/U maps. Outputs use the
`lensedcmb_new_XXXXX` prefix. Existing `lensedcmb_XXXXX.fits` files and other
input maps are preserved.

The repository includes the production, theory and validation code, together
with precomputed input spectra and installation tools. The convergence FITS
maps are external inputs held in the cluster archive. The executable code does
not depend on utility modules from `A_NEW_WORK_MAIN`.

The ten-realization acceptance results below are historical generation checks.
The manuscript's later analysis of 50 native NSIDE=2048 realizations is recorded
separately in [validation_50](../validation_50/README.md).

## Installation from GitHub

From the repository root, download the verified base runtime, seven dependency
wheels and recorded theory-source files, then install the separate CMB runtime:

```bash
python3 scripts/fetch_assets.py --cmb --theory-sources
cd cmb_sehgal2010
./setup_cmb_env.sh
./run_cmb.sh --help
```

The download helper requires a system Python 3. The installed scientific
runtime targets Linux x86-64 and Python 3.12. The numerical source code and
precomputed spectra are the same as those used in the validated production
package; the large sky maps are not downloaded by this helper.

For offline installation, transfer the release assets and theory-source files
from a networked host using the instructions in the
[root README](../README.md). Once `runtime_base.tar.gz` and `wheels/` are present
in this directory, installation requires no network access:

```bash
cd /your/final/path/cmb-sky-sim/cmb_sehgal2010
./setup_cmb_env.sh
./run_cmb.sh --help
```

The installer verifies the base archive's SHA256 and the sizes and SHA256
values of all seven wheels against `evidence/offline_wheels_manifest.json`.
It extracts the runtime into this directory's `runtime/`, runs `conda-unpack`
at the final installation prefix, installs the pinned wheels with
`pip --no-index --no-deps`, and checks dependencies and imports. The separate
Ali_sky runtime is not modified.

Keep the unpacked `runtime/` at its installation path. To use another location,
install a fresh copy there from the archived assets. Running the installer
again verifies the assets and reinstalls the same pinned wheels. An explicit
base-archive path is also accepted:

```bash
./setup_cmb_env.sh /path/to/runtime_base.tar.gz
```

The base archive's SHA256 is
`fa263b88f8f101e71d0f7de31cc68593202fa5d52ee7d5dc55b279d7dd94c577`.
The seven overrides are NumPy 1.26.4, CAMB 1.6.1, ducc0 0.38.0,
contourpy 1.3.3, SymPy 1.14.0, mpmath 1.3.0 and lenspyx 2.0.0; see
[requirements-extra.txt](requirements-extra.txt). The base runtime also
provides healpy 1.19.0, SciPy 1.14.1, Astropy and the required shared libraries.
The lenspyx wheel uses the original 2.0.0 numerical and data source files with
an adjusted packaging backend. Verification of the 43 original library files
is recorded in [the vendor manifest](evidence/lenspyx_vendor_manifest.json).

## Theory and model assumptions

[Sehgal et al. (2010), Section 1](https://arxiv.org/abs/0908.0540) gives
Ω_b = 0.044, Ω_m = 0.264, h = 0.71, n_s = 0.96, σ₈ = 0.8 and the rounded
value Ω_Λ = 0.736. This calculation uses H0 = 71, ω_b = 0.0221804 and
ω_c = 0.110902, with A_s normalized to linear σ₈(z = 0) = 0.8. The completed
calculation gives A_s = 2.1500364844398247 × 10⁻⁹ and
σ₈ = 0.8000000000000018. CAMB includes radiation when enforcing flatness,
so the calculated Ω_de is 0.73591700649.

The original authors' complete CAMB parameter file has not been recovered.
The additional assumptions used here are τ = 0.08, massless neutrinos,
N_eff = 3.046, YHe = 0.24, TCMB = 2.7255 K, a scalar pivot of 0.05/Mpc,
r = 0, zero running, Ω_K = 0, w = −1 and wa = 0. The calculation uses CAMB's
default TanhReionization/Recfast models and the original Smith2003 Halofit.
The approximate optical depth of 0.08 in the paper's gas model motivates the
chosen τ; it does not establish the exact optical depth used for the original
primary CMB. Likewise, the 2.726 K value in the authors' SED code does not
establish their CAMB background temperature. Source evidence, unrecovered
settings and the numerical configuration are documented in
[the parameter audit](theory/source_parameter_audit.md),
[the configuration](theory/sehgal2010_config.json) and
[the source manifest](theory/source_evidence/source_manifest.json).

The primary spectra in `theory/outputs/sehgal2010_theory.npz` are raw C_ell in
μK_CMB², with zero-based indices 0–6500 and columns **TT, EE, BB, TE**.
Primordial BB is exactly zero for r = 0. The primary-CMB cutoff is
ell = m = 5999, and the convergence/gradient-deflection cutoff is L = m = 6143.
The convergence `map2alm` transform uses three iterations without ring or
pixel weights. The conversion is φ_LM = 2κ_LM/[L(L+1)] for L > 0, with
φ_00 = 0 and L = 1 retained; the gradient deflection is √[L(L+1)]φ_LM.
Lenspyx remaps directly to HEALPix NSIDE=2048 with epsilon = 10⁻⁶.

Output T/Q/U maps use float64, RING ordering, the COSMO polarization convention
and μK_CMB units. Generation applies no beam, mask, pixel-window correction or
resolution downgrade. CAMB's lensing-potential spectrum is a theoretical
reference; the actual remapping uses the matched stored mmDL convergence map.

The theory outputs also contain a reference matched to the primary cutoff
5999 and potential cutoff 6143, an angular-integration sampling control and a
zero-lensing control. The full CAMB, cutoff-matched CAMB and conditional
measured-convergence references are labeled separately. A reference built
from the measured mean convergence power does not capture all higher-order
statistics of non-Gaussian convergence. Array conventions are detailed in
[the theory README](theory/README.md).

## Coordinates, seeds and output files

The producer retains the convergence map's orientation and performs no
coordinate rotation. A valid FITS `COORDSYS` header determines the output
frame. If that header is absent, the default is the dataset owner's declared
Galactic frame, G. This declaration differs from the upstream mmDL ICRS
provenance; the generator does not infer or certify an unrecorded rotation.
`--owner-coordinate C` or `E` changes only the fallback declaration for a
missing header, not the map orientation. Each JSON sidecar records the frame
provenance and this limitation.

The default master seed is 20261011. `SeedSequence` derives a uint32 seed for
each realization from a fixed domain, the master seed and the realization ID.
Scheduling order and thread count do not change that seed. This scheme does
not recover unrecorded draws from the old generator or guarantee that uint32
seeds can never collide. Set `--master-seed` explicitly to select another draw.

For ID 0, the input is `DATA_ROOT/00000/kappa_00000.fits`, and the outputs are:

- `00000/lensedcmb_new_00000.fits`
- `00000/lensedcmb_new_00000.json`
- `00000/lensedcmb_new_00000.input_spectra.npz`

Together, these files record convergence, theory and code hashes, random seeds,
run settings and full-map finite-value checks. Complete existing products can
be reused after their generation settings and hashes are verified. Partial
outputs or mismatched settings are rejected. Only explicit
`--execute --overwrite` replaces outputs with the new prefix, after preserving
them in `lensedcmb_new_backups/`. The original CMB files remain in place.

## Inspect a plan, then generate maps

The following cluster path is an example. It can be replaced by an accessible
archive with the same five-digit realization-directory structure:

```bash
CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
./run_cmb.sh --mmdl-root "$CMB_DATA_ROOT" \
  --theory "$PWD/theory/outputs/sehgal2010_theory.npz" \
  --ids 0,1 --master-seed 20261011 \
  --manifest "$PWD/evidence/example_plan.json"
```

By default, the command checks the inputs without modifying them and writes a
plan. To generate maps, use `--execute` inside an allocated compute job and set
the thread count to match the allocated CPUs:

```bash
./run_cmb.sh --mmdl-root "$CMB_DATA_ROOT" \
  --theory "$PWD/theory/outputs/sehgal2010_theory.npz" \
  --ids 0 --master-seed 20261011 --threads 16 --execute \
  --manifest "$PWD/evidence/example_production.json"
```

Realizations can also be selected with `--start-id 0 --count 10` or
`--ids-file /path/ids.txt`. The campaign-level `--manifest` must be written
outside the input data tree.

## Slurm workflow

Production runs independently for each ID and does not require MPI. This
ten-realization job template uses the selected installation path, one process
and 16 CPUs per array task, excluding aliws001, aliws003 and aliws004:

```bash
export CMB_RELEASE_ROOT=$(pwd -P)
export CMB_DATA_ROOT=/sharefs/alicpt/users/chenwz/Testarea/FORECAST/A_NEW_WORK_ILC/mmDL_maps_2048
mkdir -p logs evidence
cat > cmb_job.slurm <<'SLURM'
#!/usr/bin/env bash
#SBATCH --partition=ali
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --array=0-9%10
#SBATCH --exclude=aliws001,aliws003,aliws004
set -euo pipefail
: "${CMB_RELEASE_ROOT:?}" "${CMB_DATA_ROOT:?}"
exec "${CMB_RELEASE_ROOT}/run_cmb.sh" \
  --mmdl-root "${CMB_DATA_ROOT}" \
  --theory "${CMB_RELEASE_ROOT}/theory/outputs/sehgal2010_theory.npz" \
  --ids "${SLURM_ARRAY_TASK_ID}" --threads "${SLURM_CPUS_PER_TASK}" \
  --master-seed 20261011 --owner-coordinate G --execute \
  --manifest "${CMB_RELEASE_ROOT}/evidence/production_${SLURM_ARRAY_TASK_ID}.json"
SLURM
sbatch --export=ALL,CMB_RELEASE_ROOT,CMB_DATA_ROOT \
  --output="$CMB_RELEASE_ROOT/logs/cmb_%A_%a.log" cmb_job.slurm
```

For the remaining IDs, `sbatch --array=10-499%48` overrides the array range and
concurrency limit; retain the other arguments above. The archived scripts in
`slurm/` use the absolute installation, data and log paths of the validated
cluster runs. When relocating the package, use this template or update those
paths consistently. The supplied completed theory does not need to be
recomputed. Use a separate output directory for a new theory calculation to
preserve the release reference.

## Verification and pilot generation

After installing the environment, run the lightweight convention and
invalid-input tests:

```bash
source ./_cmb_env.sh
"$CMB_RUNTIME/bin/python" -m unittest discover -s theory -p 'test_*.py'
"$CMB_RUNTIME/bin/python" -m unittest discover -s production -p 'test_*.py'
"$CMB_RUNTIME/bin/python" -m unittest discover -s validation -p 'test_*.py'
"$CMB_RUNTIME/bin/python" validation/validate_inputs.py --root . \
  --report evidence/local_theory_input_audit.json
```

These small-array tests are separate from validation with native-resolution
production maps. The theory calculation and native NSIDE=2048 pilot for IDs
0–9 have recorded execution evidence. The independent theory review is in
[the review record](theory/source_evidence/independent_theory_review.json).
The published pilot spectra, reference statistics and binned summaries are in
[validation/pilot_summary](validation/pilot_summary/), and its plots are in
[validation/figures](validation/figures/). The full native-measurement and
summary commands are documented in [the validation design](validation/design.md).
The original cluster working directory `validation/evidence/pilot_review/`
mentioned in historical records is not included in this GitHub distribution.

The predeclared TT/EE/TE checks and the BB model comparison are reported
separately. A completed calculation does not certify the original authors'
unrecovered configuration or every scientific use of the dataset. No spectrum
amplitude was adjusted to match BB.

For the ten-realization pilot, the fitted mean amplitude offsets relative to
self-consistent CAMB are TT +0.0471%, EE −0.0525%, TE +0.2088% and BB −0.6662%.
The sub-percent BB difference remains statistically resolved. Using the mean
two-point spectrum of the actual ten convergence maps as a conditional lensing
reference gives a BB offset of −0.0395%. These fits are validation diagnostics
and were not used to modify the generated maps. See
[the quantitative results](validation/RESULTS.md) and
[the ten-realization mean and sample-standard-deviation plot](validation/figures/native_cmb_pilot_10.pdf).

The file and generation-record checks for all 500 realizations, including
preservation of the original inputs, are recorded in
[evidence/campaign_status_final.json](evidence/campaign_status_final.json).
That completeness check has a different scope from the ten-map pilot spectra
at ell < 1000. The later 50-realization manuscript comparison is available in
[validation_50](../validation_50/README.md).
