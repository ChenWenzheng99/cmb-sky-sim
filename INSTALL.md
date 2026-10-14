# Runtime installation

The release uses its own Python and MPI runtime. It does not import utility
code from the simulation archive or require an existing conda environment.
The data paths in a configuration may still point to the original cluster
maps. The bundled runtime is built for Linux x86-64; its tested host and
relocation checks are recorded in the runtime manifest.

## Install the runtime Release asset

Clone or extract the source package at its final location, download the verified
runtime with `python3 scripts/fetch_assets.py --runtime`, then run:

```bash
cd cmb-sky-sim
./setup_env.sh
```

The default archive is `runtime.tar.gz` next to this script. If the runtime
is distributed separately, give its path instead:

```bash
./setup_env.sh /path/to/ali_sky_runtime_1.0.tar.gz
```

Release 1.1.0 changes the observation source selection to the published
`lensedcmb_new_ID.fits` maps. Its dependency runtime and seven-template
cache reuse the immutable 1.0 archives. Install 1.1 at its own final
prefix; existing 1.0 source installations and archives remain separate.

The script extracts to `runtime/`, runs `conda-unpack` with the bundled
Python, checks package requirements, and imports the numerical packages and
MPI bindings. It needs neither conda nor a separate system Python. Keep the
unpacked runtime at this path; if the package is moved after this step,
extract a fresh runtime archive at the destination. This follows the
[conda-pack relocation procedure](https://conda.github.io/conda-pack/).

## Run with the matching MPI library

```bash
./run_noise.sh --config CONFIG.json
./run_observe.sh --config CONFIG.json
```

`run_noise.sh` starts the coordinator with `runtime/bin/python`; the CLI
launches its matching `runtime/bin/mpiexec`, defaulting to five ranks.
Rank-count priority is the CLI override, JSON `mpi_ranks`, then
`ALI_SKY_MPI_RANKS` with fallback five. Noise realizations are independent MPI tasks; this is
not MPI-distributed PySM smoothing. `run_observe.sh` starts the observation
dispatcher, which uses the same runtime launcher when generating noise is
required. Production calculations should run inside a suitable Slurm
allocation. The cluster submission helpers `submit_noise.sh` and
`submit_observe.sh` allocate one node and five tasks with the Ali
partition and regular QoS; their resource settings and logs are described
in README.md. The MPI runtime and `mpi4py` must come from the same installed
environment, consistent with the
[official mpi4py installation guidance](https://mpi4py.readthedocs.io/en/stable/install.html).

The scripts use only the source package on `PYTHONPATH`, disable Python
user-site packages, and default to one numerical thread per MPI rank.
`ALI_SKY_RUNTIME` can select another complete runtime directory;
`ALI_SKY_MPIEXEC` is set to that directory's MPI launcher. The wrappers set
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`NUMBA_NUM_THREADS` to one unless the caller supplied them.

## PySM templates

`PYSM_LOCAL_DATA` defaults to `template_cache/` next to the source package.
The separately fetched cache preserves the upstream relative filenames for the
fixed `d9`, `s4`, `a1` and `f1` presets. It contains seven templates;
`runtime_evidence/template_manifest.json` records their URLs, sizes and SHA256.
Fetch it with `python3 scripts/fetch_assets.py --templates`, or set
`PYSM_LOCAL_DATA` to an existing accessible provider-template mirror.
It is separate from the Python runtime and from cached frequency maps.
Set `PYSM_LOCAL_DATA` to another mirror if needed. The configured
`template_dir` may select a different template directory explicitly.

PySM loads components with `Sky` and evaluates emission at the requested
frequency; outputs are converted from its default RJ units using
thermodynamic equivalencies. See the
[official PySM interface and unit conventions](https://pysm3.readthedocs.io/en/stable/basic_use.html).
The supplied presets use fixed upstream templates; a cached map is not an
independent foreground realization.

## Rebuild the environment

`environment.yml` pins the main packages. `conda-linux-64.explicit.txt`
records the exact conda builds and dependencies used for the delivered
archive; `requirements-runtime.txt` lists the Python distributions.
Rebuild from the explicit file on a compatible Linux x86-64 host:

```bash
conda create --prefix /absolute/path/to/ali_sky_runtime --file conda-linux-64.explicit.txt
```

To resolve the pinned specification again:

```bash
conda env create --prefix /absolute/path/to/ali_sky_runtime --file environment.yml
```

Both routes include the MPI implementation. The explicit file is the
build lock; resolving the YAML later may select different dependency
builds. PySM's supported install methods and MPI extras are described in
its [official installation documentation](https://pysm3.readthedocs.io/en/),
and healpy supplies binary packages as described in its
[official installation documentation](https://healpy.readthedocs.io/en/latest/install.html).
The runtime Release asset supplies these dependencies for offline deployment
after the initial download. Template and large sky inputs remain separate.
