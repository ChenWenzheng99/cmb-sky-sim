# CMB sky simulation tools

This release generates configurable SAT/LAT noise, combines archived CMB and
foreground maps into simulated observations, and provides a Sehgal2010-based
CAMB/lenspyx CMB-generation extension. The observation/noise code is version
**1.1.0**; the CMB extension retains its validated **1.0** numerical code.

The large sky maps remain in the IHEP data archive. Input paths, coordinate
conventions, units and realization IDs are configurable. The original-code
and data license decisions remain pending; see [LICENSE_STATUS.txt](LICENSE_STATUS.txt).

## Install and run

On Linux x86-64, obtain the source repository and its checked Python/MPI runtime:

```bash
git clone https://github.com/ChenWenzheng99/cmb-sky-sim.git
cd cmb-sky-sim
python3 scripts/fetch_assets.py --runtime
./setup_env.sh
./submit_noise.sh --config examples/noise_lat.json
./submit_observe.sh --config examples/observe_lat150.json
```

The bootstrap download requires a system Python 3. The subsequent scientific
commands use the supplied Python 3.12 runtime and its matching MPI library.
The runtime is a separately downloadable GitHub Release asset, approximately
437 MB; its SHA256 is checked before installation. No existing author utility
modules or conda environment are required. Installation needs network access
only while fetching the assets. On a disconnected compute node, copy verified
assets from a networked host and use
`python3 scripts/fetch_assets.py --runtime --local-assets /path/to/assets`.
For an existing packed runtime, `./setup_env.sh /path/to/runtime.tar.gz` works.

Noise generation does not require sky maps. Observation synthesis requires the
selected archive files and their filesystem permissions. The supplied example
uses IHEP paths; edit the `repository` fields for another accessible mirror.
For PySM generation at additional frequencies, use an existing template cache
with `PYSM_LOCAL_DATA=/path/to/templates`, or fetch the recorded upstream files:

```bash
python3 scripts/fetch_assets.py --templates
```

The template download is approximately 1.26 GB. The helper preserves the
provider filenames and verifies every checksum; the template files are not
stored in GitHub. `--plan` lists the files without downloading. The cache
supports the supplied d9/s4/f1/a1 presets; other models may require additional
provider templates. Read [INSTALL.md](INSTALL.md) for runtime relocation and
MPI details, and [the CMB extension guide](cmb_sehgal2010/README.md) for remapping.

To install that separate CMB runtime and its seven dependency wheels:

```bash
python3 scripts/fetch_assets.py --cmb --theory-sources
cd cmb_sehgal2010
./setup_cmb_env.sh
./run_cmb.sh --help
```

The CMB generator needs the existing matched mmDL convergence maps. Its default
command produces a plan; `--execute` enables map generation. It does not create
additional mmDL convergence or extragalactic realizations.

## Cluster execution

On the IHEP Slurm cluster, the submission scripts request one node in
`ali` with `regular` QoS, five tasks and one CPU per task. Noise defaults
to 64 GB and 12 hours; observation defaults to 96 GB and 24 hours.
Nodes aliws001, aliws003 and aliws004 are excluded. Standard output and
errors go to `logs/` beside the source package. These are allocation
defaults, not measured memory or runtime guarantees.

Submitters forward `--config` and the same CLI arguments as the run
scripts. `ALI_SKY_JOB_RANKS`, `ALI_SKY_JOB_MEM`, `ALI_SKY_JOB_TIME`,
`ALI_SKY_JOB_PARTITION`, `ALI_SKY_JOB_QOS` and `ALI_SKY_JOB_LOG_DIR`
override allocation settings. `ALI_SKY_JOB_EXCLUDE` adds further excluded
nodes while retaining the three exclusions above. If JSON or CLI
requests more than five MPI ranks, increase the allocated task count too:

```bash
ALI_SKY_JOB_RANKS=8 ./submit_noise.sh --config examples/noise_lat.json --mpi-ranks 8
```

`ALI_SKY_JOB_RANKS` controls Slurm `--ntasks`, while JSON/CLI `mpi_ranks`
controls the noise launcher. Keep them consistent. Inside an existing
compute allocation, call the run scripts directly:

```bash
./run_noise.sh --config examples/noise_lat.json
./run_observe.sh --config examples/observe_lat150.json
```

Both launchers start a serial coordinator. The noise CLI launches its
matching MPI workers; do not wrap the observation coordinator in
`mpiexec`. Rank-count priority is `--mpi-ranks`, JSON `mpi_ranks`, then
`ALI_SKY_MPI_RANKS` with fallback five. A noise campaign assigns independent
(frequency, realization) tasks to MPI ranks. Observation synthesis follows
the noise stage serially. Allocate enough memory for native NSIDE2048
maps and the requested concurrent workers; the numerical thread default
is one per worker.

To inspect a configuration without producing maps:

```bash
source ./_runtime_env.sh
"$ALI_SKY_RUNTIME/bin/python" -m ali_sky.cli validate-config --config examples/observe_lat150.json --mode observe
"$ALI_SKY_RUNTIME/bin/python" -m ali_sky.cli validate-config --config examples/noise_lat.json --mode noise
```

Noise and observation commands accept `--n-realizations`, `--start-id`,
`--seed`, `--output-dir` and `--mpi-ranks`. CLI values override JSON values.
Changing count or start ID replaces an explicit JSON `ids` list.
Paths in JSON, including relative CLI output paths, resolve against the
JSON file's directory. The examples place products under `outputs/` next
to this README.

## Example configurations

| Configuration | Action |
|---|---|
| [noise_sat.json](examples/noise_sat.json) | Generate 50 NSIDE2048 SAT realizations in each of six bands. |
| [noise_lat.json](examples/noise_lat.json) | Generate 50 NSIDE2048 LAT realizations in each of six bands. |
| [observe_lat150.json](examples/observe_lat150.json) | Combine archive sky ID 0 and archive LAT150 noise ID 0 with the default mask. |
| [custom_95_150.json](examples/custom_95_150.json) | Generate two custom noise bands automatically, evaluate/cache Galactic emission at 95 GHz, and combine sky ID 0 at 95 and 150 GHz. |
| [with_mask.json](examples/with_mask.json) | Select an explicit native mask path and its C input frame. |
| [use_existing_noise.json](examples/use_existing_noise.json) | Read the FITS noise products from `noise_lat.json` instead of generating or reading the legacy NPY maps. |

The custom example's depths and beams are illustrative user inputs,
including white-only noise. They are not a calibrated 95 GHz instrument
forecast. Set the parameters appropriate to the intended instrument.
Run `noise_lat.json` before `use_existing_noise.json`, or change its
`noise_dir` and pattern to a prepared noise collection.

The default archive root is
`/sharefs/alicpt/users/chenwz/Testarea/FORECAST`. On another installation,
set `repository.sky_dir`, `repository.galactic_dir`,
`repository.noise_dir` and `repository.mask` to that installation's files.
Only data files are read from these paths. Sky IDs must exist explicitly;
the CMB reader and preflight require
`SKY_DIR/ID/lensedcmb_new_ID.fits` (five-digit IDs). These are the 500
Sehgal2010-CAMB-primary maps remapped with the corresponding mmDL kappa,
with zero primordial tensors. Retained `lensedcmb_ID.fits` maps are ignored;
there is no legacy CMB fallback. Each selected source path and SHA256 is
recorded in the observation sidecar. The independent CMB-generation
extension supplies its own theory, seeds, signatures and installation
instructions; observation synthesis reads its published FITS products.
Noise IDs are matched explicitly through `noise_ids` and never wrapped
modulo 50. The archive's existing 50 noise maps are separate from the
seeded products generated by this package: their original random seeds
are unavailable, so regenerating their exact pixels is not claimed.

## Coordinates, maps and masks

The default output frame is Galactic (`coordinate: "G"`). Following the
dataset owner's declaration on 2026-10-11, the default missing-header
assumptions for headerless XFG and legacy noise are G. The published
`lensedcmb_new_ID.fits` CMB products explicitly carry G headers. This declaration
does not independently verify those files' coordinates; historical
provider/source-trace evidence identifies upstream mmDL products as
C/ICRS and remains a separate provenance record. Input FITS `COORDSYS`
headers take precedence over a missing-header assumption.

The component-separated Galactic maps are G. With G inputs and G output,
no map rotation is performed. For a collection known to contain C maps,
set `repository.sky_coordinate: "C"`; independently set the appropriate
noise frame. Selecting output C rotates all G inputs, including Galactic
maps, to C. Map transfer uses native-grid interpolation and the spin-2
tangent-basis rotation for Q/U. It is not a spherical-harmonic rotation
or a resolution downgrade. Scalar extragalactic intensity follows the
same input/output frame choice.

The default archive mask is explicitly C:
`A_NEW_WORK_MAIN/sim/mask_c/AliCPT_LAT_cross_PLANCK_2048_bi_C.fits`.
It is transferred to G by nearest-pixel lookup for the default output.
For another mask, edit `mask` and declare `mask_coordinate` as G or C
as appropriate. A FITS coordinate header still takes precedence.
Mask weights must lie in [0,1], and input resolution must equal `nside`.
Set `use_repository_mask: false` and omit `mask` for full-sky output.
The package applies the selected mask directly without apodization.

All generated maps are float64 HEALPix RING T/Q/U, with thermodynamic
`uK_CMB`, explicit `COORDSYS` and `POLCCONV=COSMO`, FITS checksums, and a
JSON sidecar. The readers check resolution, finite values, HEALPix UNSEEN,
units and available frame metadata. No automatic input downgrade occurs.

## Noise model and reproducibility

The default bands are 27, 40, 90, 150, 225 and 280 GHz. Their SAT/LAT
depths, beams and red-noise parameters are defined in
[ali_sky/noise.py](ali_sky/noise.py). An unlisted frequency requires a
`channels` entry with `frequency_GHz`, `nlev_uK_arcmin` and
`beam_fwhm_arcmin`; red-noise parameters can be specified explicitly.
`include_red: false` selects white-only noise.

| Channel parameter | Meaning and units |
|---|---|
| `nlev_uK_arcmin` | Temperature white-noise depth, thermodynamic microkelvin arcmin; scalar Q/U white C_ell is twice the T value. |
| `N_red_T` | Temperature red C_ell amplitude in thermodynamic microkelvin squared. |
| `n_red_P` | Polarization red-noise amplitude in microkelvin arcmin; null uses the white depth. |
| `l_knee_T`, `l_knee_P` | Positive multipole scales for the red terms. |
| `alpha_knee_T`, `alpha_knee_P` | Finite power-law exponents. |
| `beam_fwhm_arcmin` | Gaussian sky beam FWHM in arcmin, recorded but not applied by the noise stage. |

T, Q and U are synthesized as three independent scalar fields, preserving
the archived generator convention rather than treating Q/U as E/B inputs.
The model's array entry at harmonic index ell is evaluated at ell+1;
`scale_factor: 1` preserves the original same-grid generation convention.
The noise stage applies neither the sky beam nor a mask or pixel window.
New realization seeds depend on the base seed, instrument, frequency and
ID, so task order and MPI rank do not choose different random streams.
Existing matching generated products may be reused only after their
sidecar configuration and hash checks; observations refuse an existing
output unless `overwrite` is explicitly true.

Noise products are
`OUTPUT/INSTRUMENT/FREQUENCY/noise_00000.fits` and `.json`, with
`noise_configuration.json` and `noise_campaign.json` at the campaign root.
The reader also supports legacy `TQU_N_ID.npy` arrays and verifies their
complete RING index column. FITS noise can be selected with `noise_dir`,
`noise_pattern` and `noise_coordinate`. Pattern fields are `instrument`,
`frequency`, `frequency_GHz` and `id`.

## Forward observations and arbitrary frequencies

The observation is
`mask * [beam * (CMB + Galactic + extragalactic) + noise]`.
The Gaussian beam acts on the combined sky's T/E/B harmonics, with
default `lmax=3*nside-1` and `harmonic_iter=1`, and is not applied to
noise. No pixel-window factor or monopole/dipole removal is applied.
The output is a map; the package does not estimate or decouple its masked
power spectrum.

For a custom transfer function, a channel's `beam_transfer_path` selects
a finite NPY vector B_ell shared by T/E/B, or a 3-by-(L+1) array ordered
T/E/B, covering the configured `lmax`. This replaces the Gaussian beam
operation; the channel still requires `beam_fwhm_arcmin` for its schema.

Galactic presets default to dust d9, synchrotron s4, free-free f1 and AME
a1. A matching component file is read from `repository.galactic_dir`;
otherwise PySM evaluates it at the requested monochromatic frequency and
caches a G-frame `uK_CMB` map. Cache identity includes the component,
preset, frequency, NSIDE and PySM version, with a file hash and template
provenance. The bundled seven templates support these fixed presets at
NSIDE2048. They are not independent foreground draws; other presets or
resolutions may require additional upstream templates.

For tSZ, CIB and radio, a genuine archived nominal-frequency map is used
when present. Otherwise the original 148 GHz map is scaled using the
retained production monochromatic SED in [ali_sky/sed.py](ali_sky/sed.py),
including the signed tSZ response. CMB and kSZ have no frequency scaling
in thermodynamic units. The four XFG components contribute T only;
extragalactic Q/U are zero. The lensing convergence map is not an additive
observation component.

With custom noise parameters and no explicit `noise_dir`, `observe`
generates all requested frequency/noise-ID combinations using its MPI
stage in `OUTPUT/generated_noise`, then combines the maps. It resolves
required sky and mask inputs before launching this stage. With explicit
`noise_dir`, it reads those maps and checks their existence instead.

Observation products are
`OUTPUT/FREQUENCY/observation_00000_noise00000_G.fits` and `.json`, plus
`observation_campaign.json`. Sidecars record source hashes, source and
output frames, beam, mask, SED factors, presets and execution versions.

## Scoped checks

Current code-version checks are recorded in `release_validation.json` and
`provenance/release_1.1_source_validation.json`. The bundled 1.0 runtime
manifest and historical native-integration reports retain their original
tested source hashes; they do not certify the 1.1 CMB filename selection.
The dependency and relocation checks are recorded in
`runtime_manifest.json`. The source tests include independently generated
small maps and noise statistics; they are not downgraded NSIDE2048
archive maps. To run the complete source test suite with the delivered runtime:

```bash
source ./_runtime_env.sh
"$ALI_SKY_RUNTIME/bin/python" -m unittest discover -s tests
```

These checks test the stated software contracts. They do not reproduce
the unrecoverable legacy random seeds or certify every scientific use of
the input archive.
