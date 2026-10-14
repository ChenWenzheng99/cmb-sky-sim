#!/usr/bin/env python3
"""Native-grid Sehgal2010-primary CMB remapping using existing mmDL kappa.

The retained production cutoffs are primary lmax=mmax=5999 and kappa /
deflection lmax=mmax=6143 at NSIDE=2048. No coordinate rotation, beam,
mask, pixel-window correction, or map downgrade is performed. The theory
archive supplies raw thermodynamic C_l in explicit TT, EE, BB, TE order;
its CAMB phi spectrum is a reference and is NOT used as the deflection.

The command defaults to a read-only plan. Map writing requires --execute.
The original producer is preserved in original_snapshot/ with its hash.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import gc
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
from pathlib import Path
import platform
import socket
import sys
import tempfile
import threading
import time
from typing import Any


SCHEMA = "sehgal2010_mmdl_native_cmb_v1"
THEORY_SCHEMA = "sehgal2010_camb_theory_v1"
NSIDE = 2048
PRIMARY_LMAX = 5999
DEFLECTION_LMAX = 6143
KAPPA_ITER = 3
EPSILON = 1e-6
COLUMNS = ["TT", "EE", "BB", "TE"]
ORIGINAL_SHA256 = "ab7c9bb4dfc0f1bcdfbbb04cf7658a8feb73837193618a8013067111de141d21"
SEED_POLICY = "numpy.SeedSequence(domain_sha256_u32, master_seed_u64, realization_id_u64); uint32_v1"
SEED_DOMAIN = "sehgal2010-primary-mmdl-kappa-native2048-v1"
_RNG_LOCK = threading.Lock()


def sha256(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def file_record(path: Path | str) -> dict[str, Any]:
    path = Path(path).resolve(strict=True)
    before = path.stat()
    digest = sha256(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f"File changed while hashing: {path}")
    return {"path": str(path), "bytes": after.st_size,
            "mtime_ns": after.st_mtime_ns, "sha256": digest}


def checked_integer(value: Any, name: str, minimum: int = 0,
                    maximum: int = 2 ** 64 - 1) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def realization_seed(master_seed: int, realization_id: int) -> int:
    """Stable task seed, independent of scheduling, list order and rank count.

    This is a new reproducible policy. It does not recover the original
    producer's unrecorded seeds, and uint32 seeds are not collision-proof.
    """
    import numpy as np
    checked_integer(master_seed, "master_seed")
    checked_integer(realization_id, "realization_id")
    domain = hashlib.sha256(SEED_DOMAIN.encode()).digest()
    entropy = [int.from_bytes(domain[i:i + 4], "little") for i in range(0, 32, 4)]
    entropy += [master_seed & 0xffffffff, master_seed >> 32,
                realization_id & 0xffffffff, realization_id >> 32]
    return int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint32)[0])


def full_statistics(array: Any, *, expected_shape: tuple[int, ...] | None = None) -> dict[str, Any]:
    """Full-array finite and UNSEEN checks without a second full-size copy."""
    import numpy as np
    values = np.asarray(array)
    if expected_shape is not None and values.shape != expected_shape:
        raise ValueError(f"Expected shape {expected_shape}; got {values.shape}")
    if values.dtype.kind not in "fi":
        raise ValueError("Expected real numeric maps")
    flat = values.reshape(-1)
    finite_count = unseen_count = 0
    minimum, maximum = math.inf, -math.inf
    total = total_sq = 0.0
    # healpy.UNSEEN=-1.6375e30. isclose follows healpy.mask_bad's relative tolerance.
    for start in range(0, flat.size, 1024 * 1024):
        chunk = flat[start:start + 1024 * 1024]
        finite = np.isfinite(chunk)
        finite_count += int(np.count_nonzero(finite))
        unseen_count += int(np.count_nonzero(np.isclose(chunk, -1.6375e30,
                                                       rtol=1e-5, atol=1e-8)))
        if not bool(np.all(finite)):
            raise ValueError("Non-finite values in full native array")
        minimum = min(minimum, float(np.min(chunk)))
        maximum = max(maximum, float(np.max(chunk)))
        total += float(np.sum(chunk, dtype=np.float64))
        total_sq += float(np.dot(chunk.astype(np.float64, copy=False),
                                 chunk.astype(np.float64, copy=False)))
    if unseen_count:
        raise ValueError(f"UNSEEN values in full native array: {unseen_count}")
    if flat.size == 0:
        raise ValueError("Empty array")
    mean = total / flat.size
    return {"shape": list(values.shape), "dtype": values.dtype.name, "count": int(flat.size),
            "finite_count": finite_count, "unseen_count": unseen_count,
            "min": minimum, "max": maximum, "mean": mean,
            "population_std": math.sqrt(max(0.0, total_sq / flat.size - mean * mean))}


def normalize_coordinate(value: str) -> str:
    normalized = str(value).strip().upper()
    choices = {"G": "G", "GALACTIC": "G", "C": "C", "EQUATORIAL": "C",
               "CELESTIAL": "C", "ICRS": "C", "E": "E", "ECLIPTIC": "E"}
    if normalized not in choices:
        raise ValueError(f"Unrecognized coordinate frame: {value!r}")
    return choices[normalized]


def coordinate_provenance(header: dict[str, Any], owner_coordinate: str = "G") -> dict[str, Any]:
    declared = normalize_coordinate(owner_coordinate)
    value = header.get("COORDSYS")
    if value is not None and str(value).strip():
        chosen = normalize_coordinate(value)
        origin = "source FITS COORDSYS header"
        verified = True
    else:
        chosen = declared
        origin = "owner-declared frame for headerless repository kappa"
        verified = False
    return {"output_coordinate": chosen, "source_header_COORDSYS": value,
            "owner_declared_coordinate": declared, "selection": origin,
            "source_frame_available_in_header": verified,
            "coordinate_rotation_performed": False,
            "directions": "Existing kappa directions retained; the generated primary is remapped in the same frame",
            "caveat": ("A headerless local kappa frame is not physically re-certified by this producer. "
                       "The owner's G declaration differs from the upstream native mmDL ICRS provenance; "
                       "no unrecorded intermediate rotation is inferred.") if not verified else
                      "Source header takes priority; directions are retained without rotation"}


def read_kappa_header(path: Path | str) -> dict[str, Any]:
    from astropy.io import fits
    with fits.open(path, memmap=True, checksum=True) as hdus:
        header = hdus[1].header
        if header.get("NSIDE") != NSIDE or str(header.get("ORDERING", "")).upper() != "RING":
            raise ValueError("Kappa must be source NSIDE2048 RING; no conversion is allowed")
        if header.get("TFIELDS") != 1:
            raise ValueError("Expected one scalar kappa field")
        if header.get("INDXSCHM", "IMPLICIT").upper() != "IMPLICIT":
            raise ValueError("Expected an implicit full-sky kappa FITS table")
        if header.get("FIRSTPIX", 0) != 0 or header.get("LASTPIX", 12 * NSIDE ** 2 - 1) != 12 * NSIDE ** 2 - 1:
            raise ValueError("Kappa pixel range must cover the full native sky")
        unit = header.get("TUNIT1")
        if unit is not None and str(unit).strip().lower() not in ["", "1", "dimensionless"]:
            raise ValueError(f"Kappa must be dimensionless, not {unit!r}")
        return {name: header.get(name) for name in ["NSIDE", "ORDERING", "TFIELDS", "INDXSCHM",
                "FIRSTPIX", "LASTPIX", "COORDSYS", "TUNIT1", "POLCCONV", "TTYPE1"]}


def load_theory(archive_path: Path | str, manifest_path: Path | str | None = None) -> dict[str, Any]:
    """Read a hash-certified raw-C_l archive; never guess column conventions."""
    import numpy as np
    archive_path = Path(archive_path).resolve(strict=True)
    manifest_path = Path(manifest_path or archive_path.with_name("theory_manifest.json")).resolve(strict=True)
    archive_record, manifest_record = file_record(archive_path), file_record(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != THEORY_SCHEMA or manifest.get("status") != "COMPLETE_THEORY_WITH_DECLARED_ASSUMPTIONS":
        raise ValueError("Expected a completed Sehgal2010 theory manifest")
    expected = manifest.get("artifact_manifest", {}).get(archive_path.name)
    if expected is None or expected.get("sha256") != archive_record["sha256"] or expected.get("bytes") != archive_record["bytes"]:
        raise ValueError("Theory archive hash/size does not match its manifest")
    cutoffs = manifest.get("producer_cutoffs", {})
    if any(cutoffs.get(k) != v for k, v in {"primary_lmax": PRIMARY_LMAX,
           "kappa_deflection_lmax": DEFLECTION_LMAX, "nside": NSIDE}.items()):
        raise ValueError("Theory manifest producer cutoffs differ from the retained source settings")
    contract = manifest.get("array_contract", {})
    arrays: dict[str, Any] = {}
    with np.load(archive_path, allow_pickle=False) as archive:
        for key in ["ell", "unlensed_raw_cl", "lensed_raw_cl", "producer_unlensed_raw_cl", "phi_raw_cl", "phi_PP"]:
            if key not in archive:
                raise ValueError(f"Missing theory array {key}")
            arrays[key] = np.array(archive[key], copy=True)
    ell = arrays["ell"]
    if ell.ndim != 1 or ell.dtype != np.dtype("int64") or len(ell) <= DEFLECTION_LMAX or not np.array_equal(ell, np.arange(len(ell))):
        raise ValueError("Theory ell must explicitly cover contiguous multipoles from zero through 6143")
    for key, array in arrays.items():
        if not np.all(np.isfinite(array)):
            raise ValueError(f"Non-finite theory spectrum: {key}")
        description = contract.get(key, {})
        if description.get("shape") != list(array.shape) or description.get("dtype") != str(array.dtype):
            raise ValueError(f"Theory array shape/dtype differs from manifest: {key}")
    for key in ["unlensed_raw_cl", "lensed_raw_cl", "producer_unlensed_raw_cl"]:
        array = arrays[key]
        if array.shape != (len(ell), 4) or array.dtype != np.dtype("float64") or contract[key].get("columns") != COLUMNS or contract[key].get("units") != "uK_CMB^2":
            raise ValueError(f"Expected TT/EE/BB/TE raw uK_CMB^2 spectrum: {key}")
        if np.any(array[:, :3] < 0):
            raise ValueError(f"Negative theory auto-spectrum: {key}")
        tolerance = 1e-12 * np.maximum(array[:, 0] * array[:, 1], np.finfo(float).tiny)
        if np.any(array[:, 3] ** 2 > array[:, 0] * array[:, 1] + tolerance):
            raise ValueError(f"Theory TT/EE/TE covariance is not positive semidefinite: {key}")
    if not np.array_equal(arrays["unlensed_raw_cl"][:, 2], np.zeros(len(ell))):
        raise ValueError("Scalar-only unlensed BB must be exactly zero")
    primary = arrays["producer_unlensed_raw_cl"]
    if not np.array_equal(primary[:PRIMARY_LMAX + 1], arrays["unlensed_raw_cl"][:PRIMARY_LMAX + 1]) or np.any(primary[PRIMARY_LMAX + 1:] != 0):
        raise ValueError("Producer primary must retain raw C_l through 5999 and be exactly zero above it")
    for key in ["phi_raw_cl", "phi_PP"]:
        if arrays[key].shape != (len(ell),) or arrays[key].dtype != np.dtype("float64") or contract[key].get("units") != "dimensionless" or np.any(arrays[key] < 0):
            raise ValueError(f"Expected nonnegative dimensionless phi reference: {key}")
    factor = (ell.astype(float) * (ell + 1)) ** 2 / (2 * np.pi)
    if not np.allclose(arrays["phi_PP"], factor * arrays["phi_raw_cl"], rtol=1e-12, atol=1e-30):
        raise ValueError("Theoretical phi_PP normalization mismatch")
    # These small arrays are retained; no theoretical phi is sent to remapping.
    return {"arrays": arrays, "archive": archive_record, "manifest": manifest_record,
            "configuration": manifest.get("configuration"),
            "configuration_canonical_sha256": manifest.get("configuration_canonical_sha256"),
            "normalization": manifest.get("normalization"), "versions": manifest.get("versions")}


def harmonic_transfers(lmax: int = DEFLECTION_LMAX) -> tuple[Any, Any]:
    import numpy as np
    checked_integer(lmax, "lmax", minimum=1)
    ell = np.arange(lmax + 1, dtype=np.float64)
    phi = np.zeros(lmax + 1, dtype=np.float64)
    phi[1:] = 2.0 / (ell[1:] * (ell[1:] + 1))
    deflection = phi * np.sqrt(ell * (ell + 1))
    return phi, deflection


def draw_unlensed_alms(raw_cl: Any, seed: int, *, lmax: int = PRIMARY_LMAX) -> Any:
    """Draw TT/EE/BB/TE with explicit healpy new=True, preserving global RNG."""
    import healpy as hp
    import numpy as np
    checked_integer(seed, "seed", maximum=2 ** 32 - 1)
    raw = np.asarray(raw_cl, dtype=np.float64)
    if raw.ndim != 2 or raw.shape[0] < lmax + 1 or raw.shape[1] != 4:
        raise ValueError("Expected raw C_l columns TT, EE, BB, TE")
    with _RNG_LOCK:
        previous = np.random.get_state()
        try:
            np.random.seed(seed)
            alms = hp.synalm([raw[:lmax + 1, i] for i in range(4)], lmax=lmax,
                            mmax=lmax, new=True)
        finally:
            np.random.set_state(previous)
    return alms


def input_spectrum_arrays(primary_alms: Any, kappa_alm: Any,
                          phi_alm: Any, deflection_alm: Any, raw_cl: Any) -> dict[str, Any]:
    """Existing alms only; this audit adds no map-to-alm transforms."""
    import healpy as hp
    import numpy as np
    sampled = np.asarray([hp.alm2cl(primary_alms[0]), hp.alm2cl(primary_alms[1]),
                          hp.alm2cl(primary_alms[2]), hp.alm2cl(primary_alms[0], primary_alms[1])])
    phi_transfer, deflection_transfer = harmonic_transfers()
    phi_cl = hp.alm2cl(phi_alm)
    ell_deflection = np.arange(DEFLECTION_LMAX + 1, dtype=np.int64)
    return {"ell_primary": np.arange(PRIMARY_LMAX + 1, dtype=np.int64),
            "unlensed_input_C_ell": np.asarray(raw_cl[:PRIMARY_LMAX + 1], dtype=np.float64).T.copy(),
            "unlensed_sampled_C_ell": sampled,
            "primary_columns": np.asarray(COLUMNS),
            "ell_kappa_phi_deflection": ell_deflection,
            "kappa_C_ell": hp.alm2cl(kappa_alm), "phi_C_ell": phi_cl,
            "phi_PP": phi_cl * (ell_deflection * (ell_deflection + 1)) ** 2 / (2 * np.pi),
            "deflection_C_ell": hp.alm2cl(deflection_alm),
            "kappa_phi_C_ell": hp.alm2cl(kappa_alm, phi_alm),
            "kappa_deflection_C_ell": hp.alm2cl(kappa_alm, deflection_alm),
            "phi_deflection_C_ell": hp.alm2cl(phi_alm, deflection_alm),
            "phi_transfer": phi_transfer, "deflection_transfer": deflection_transfer}


def output_paths(root: Path | str, realization_id: int) -> dict[str, Path]:
    checked_integer(realization_id, "realization_id", maximum=99999)
    directory = Path(root).resolve() / f"{realization_id:05d}"
    stem = f"lensedcmb_new_{realization_id:05d}"
    return {"directory": directory, "fits": directory / (stem + ".fits"),
            "json": directory / (stem + ".json"),
            "spectra": directory / (stem + ".input_spectra.npz"),
            "lock": directory / (stem + ".lock"),
            "kappa": directory / f"kappa_{realization_id:05d}.fits",
            "old_cmb": directory / f"lensedcmb_{realization_id:05d}.fits"}


def generation_parameters() -> dict[str, Any]:
    return {"nside": NSIDE, "ordering": "RING", "primary_lmax": PRIMARY_LMAX,
            "primary_mmax": PRIMARY_LMAX, "primary_columns": COLUMNS,
            "primary_units": "uK_CMB^2", "synalm_new": True,
            "kappa_lmax": DEFLECTION_LMAX, "kappa_mmax": DEFLECTION_LMAX,
            "kappa_map2alm_iter": KAPPA_ITER, "kappa_map2alm_pol": False,
            "kappa_use_weights": False, "kappa_use_pixel_weights": False,
            "phi_transfer": "2/[L(L+1)] for L>=1; zero at L=0; L=1 retained",
            "gradient_deflection_transfer": "sqrt[L(L+1)] times phi",
            "deflection_lmax": DEFLECTION_LMAX, "epsilon": EPSILON,
            "geometry": ["healpix", {"nside": NSIDE}], "remapping_pol": True,
            "maps_dtype": "float64", "maps_units": "uK_CMB", "polarization_convention": "COSMO",
            "beam_applied": False, "mask_applied": False, "pixel_window_applied": False,
            "downgrade_performed": False, "coordinate_rotation_performed": False}


def production_source_records() -> dict[str, Any]:
    directory = Path(__file__).resolve().parent
    snapshot = directory / "original_snapshot" / "lensed_CMB_generation.py"
    manifest = directory / "original_snapshot" / "source_manifest.json"
    record = file_record(snapshot)
    if record["sha256"] != ORIGINAL_SHA256:
        raise ValueError("Immutable original producer snapshot hash mismatch")
    original_manifest = json.loads(manifest.read_text())
    if original_manifest.get("snapshot_sha256") != ORIGINAL_SHA256:
        raise ValueError("Original source manifest mismatch")
    return {"original_snapshot": record, "original_manifest": file_record(manifest),
            "original_source": original_manifest, "release_producer": file_record(__file__)}


def make_contract(realization_id: int, master_seed: int, source: dict[str, Any],
                  theory: dict[str, Any], coordinate: dict[str, Any]) -> dict[str, Any]:
    # Runtime paths, file timestamps and rank allocation are deliberately not
    # seed inputs. Content hashes define scientific input identity for resume.
    return {"schema": SCHEMA, "realization_id": realization_id,
            "master_seed": master_seed, "seed": realization_seed(master_seed, realization_id),
            "seed_policy": SEED_POLICY, "source_kappa_sha256": source["sha256"],
            "source_kappa_bytes": source["bytes"], "theory_archive_sha256": theory["archive"]["sha256"],
            "theory_manifest_sha256": theory["manifest"]["sha256"],
            "parameters": generation_parameters(), "coordinate_provenance": coordinate,
            "runtime_thread_count_affects_seed": False,
            "producer_sha256": sha256(__file__), "original_producer_sha256": ORIGINAL_SHA256}


def _read_output_for_verification(path: Path, coordinate: str) -> list[dict[str, Any]]:
    import numpy as np
    from astropy.io import fits
    with fits.open(path, memmap=True, checksum=True) as hdus:
        header = hdus[1].header
        expected = {"NSIDE": NSIDE, "ORDERING": "RING", "COORDSYS": coordinate,
                    "POLCCONV": "COSMO", "TFIELDS": 3, "FIRSTPIX": 0,
                    "LASTPIX": 12 * NSIDE ** 2 - 1, "INDXSCHM": "IMPLICIT"}
        if any(header.get(k) != v for k, v in expected.items()):
            raise FileExistsError("Existing output FITS header differs from contract")
        if any(header.get(f"TUNIT{i}") != "uK_CMB" or not str(header.get(f"TFORM{i}", "")).endswith("D") for i in [1, 2, 3]):
            raise FileExistsError("Existing output must be float64 TQU in uK_CMB")
        if [header.get(f"TTYPE{i}") for i in [1, 2, 3]] != ["T", "Q", "U"]:
            raise FileExistsError("Existing output does not have standard TQU columns")
        for hdu in hdus:
            if "CHECKSUM" not in hdu.header or hdu.verify_checksum() != 1 or hdu.verify_datasum() != 1:
                raise FileExistsError("Existing output FITS checksum is missing or invalid")
        statistics = [full_statistics(np.asarray(hdus[1].data.field(i)).reshape(-1),
                                      expected_shape=(12 * NSIDE ** 2,)) for i in range(3)]
    return statistics


def verify_existing(paths: dict[str, Path], contract: dict[str, Any]) -> dict[str, Any] | None:
    present = [paths[key].exists() for key in ["fits", "spectra", "json"]]
    if not any(present):
        return None
    if not all(present) or any(paths[key].is_symlink() for key in ["fits", "spectra", "json"]):
        raise FileExistsError("Partial or symbolic-link output exists; explicit backup overwrite is required")
    try:
        metadata = json.loads(paths["json"].read_text())
    except (OSError, ValueError) as error:
        raise FileExistsError("Existing output sidecar is unreadable or malformed") from error
    if metadata.get("status") != "COMPLETE" or metadata.get("contract") != contract or metadata.get("contract_sha256") != canonical_sha256(contract):
        raise FileExistsError("Existing output provenance does not match requested inputs/parameters/seeds")
    for key, name in [("fits", "output"), ("spectra", "input_spectra")]:
        current = file_record(paths[key])
        expected = metadata.get(name, {})
        if current["sha256"] != expected.get("sha256") or current["bytes"] != expected.get("bytes"):
            raise FileExistsError(f"Existing output content hash mismatch: {paths[key]}")
    try:
        statistics = _read_output_for_verification(paths["fits"], contract["coordinate_provenance"]["output_coordinate"])
    except (OSError, ValueError) as error:
        raise FileExistsError("Existing output FITS failed verification") from error
    # Saved statistics are for the actual float64 serialized arrays, and must
    # agree on every full-array field before resuming.
    if statistics != metadata.get("finite_statistics", {}).get("output_serialized_TQU"):
        raise FileExistsError("Existing output full-array statistics differ from sidecar")
    return {**metadata, "reused": True, "resume_checked_at_utc": utc_now()}


@contextmanager
def exclusive_task_lock(path: Path):
    # Never remove a pre-existing lock: it may belong to an active worker.
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump({"host": socket.gethostname(), "pid": os.getpid(), "created_utc": utc_now()}, handle)
            handle.flush()
            os.fsync(handle.fileno())
        yield
    finally:
        path.unlink(missing_ok=True)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_publish(staged: Path, target: Path, overwrite: bool = False) -> None:
    """Publish on the same filesystem without a check-then-replace race."""
    if overwrite:
        os.replace(staged, target)
    else:
        os.link(staged, target)  # fails atomically if target already exists
        staged.unlink()


def backup_existing(paths: dict[str, Path]) -> dict[str, Any] | None:
    existing = [paths[k] for k in ["fits", "spectra", "json"] if paths[k].exists()]
    if not existing:
        return None
    if any(p.is_symlink() or not p.is_file() for p in existing):
        raise FileExistsError("Refusing to overwrite symbolic links or non-regular files")
    backup_dir = paths["directory"] / "lensedcmb_new_backups" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    backup_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for path in existing:
        backup = backup_dir / path.name
        os.link(path, backup)
        records.append({"original": str(path), "backup": file_record(backup)})
    report = {"created_utc": utc_now(), "strategy": "hard links preserving old inode before atomic replace",
              "files": records, "old_lensedcmb_file_touched": False}
    (backup_dir / "backup_manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def installed_versions(lenspyx: Any) -> dict[str, Any]:
    import healpy
    import numpy
    import astropy
    result = {"python": platform.python_version(), "numpy": numpy.__version__,
              "healpy": healpy.__version__, "astropy": astropy.__version__,
              "lenspyx": getattr(lenspyx, "__version__", None)}
    for package in ["lenspyx", "ducc0", "camb", "mpi4py"]:
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result.setdefault(package, "not installed")
    result["lenspyx_module_path"] = str(Path(lenspyx.__file__).resolve())
    result["lenspyx_lensing_source"] = file_record(inspect.getsourcefile(lenspyx.alm2lenmap))
    return result


def configure_threads(threads: int) -> None:
    checked_integer(threads, "threads", minimum=1, maximum=4096)
    existing = os.environ.get("OMP_NUM_THREADS")
    if "healpy" in sys.modules and existing != str(threads):
        raise RuntimeError("Set OMP_NUM_THREADS before importing healpy; it must match --threads")
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")


def generate_one(root: Path | str, theory: dict[str, Any], realization_id: int,
                 master_seed: int, *, threads: int = 16, owner_coordinate: str = "G",
                 overwrite: bool = False) -> dict[str, Any]:
    """Generate one complete native map, or verify/reuse its exact counterpart."""
    configure_threads(threads)
    import healpy as hp
    import lenspyx
    import numpy as np
    from astropy.io import fits

    started = time.monotonic()
    paths = output_paths(root, realization_id)
    if not paths["directory"].is_dir():
        raise FileNotFoundError("An existing mmDL realization directory is required")
    source = file_record(paths["kappa"])
    header = read_kappa_header(paths["kappa"])
    coordinate = coordinate_provenance(header, owner_coordinate)
    contract = make_contract(realization_id, master_seed, source, theory, coordinate)
    source_records = production_source_records()
    with exclusive_task_lock(paths["lock"]):
        if not overwrite:
            reused = verify_existing(paths, contract)
            if reused is not None:
                return reused
        # Do not back up or mutate existing files until the new output has
        # passed generation, serialization and full-source identity checks.
        kappa = hp.read_map(paths["kappa"], field=0, dtype=np.float64, nest=False)
        source_statistics = full_statistics(kappa, expected_shape=(12 * NSIDE ** 2,))
        if file_record(paths["kappa"]) != source:
            raise RuntimeError("Kappa source changed while reading")
        primary = theory["arrays"]["producer_unlensed_raw_cl"]
        primary_alms = draw_unlensed_alms(primary, contract["seed"])
        kappa_alm = hp.map2alm(kappa, lmax=DEFLECTION_LMAX, mmax=DEFLECTION_LMAX,
                              iter=KAPPA_ITER, pol=False, use_weights=False,
                              use_pixel_weights=False)
        del kappa
        gc.collect()
        phi_transfer, deflection_transfer = harmonic_transfers()
        phi_alm = hp.almxfl(kappa_alm, phi_transfer, inplace=False)
        deflection_alm = hp.almxfl(kappa_alm, deflection_transfer, inplace=False)
        spectrum_arrays = input_spectrum_arrays(primary_alms, kappa_alm, phi_alm,
                                                deflection_alm, primary)
        spectrum_arrays.update({"seed": np.asarray(contract["seed"], dtype=np.uint32),
                                "realization_id": np.asarray(realization_id, dtype=np.int64),
                                "source_kappa_sha256": np.asarray(source["sha256"]),
                                "theory_archive_sha256": np.asarray(theory["archive"]["sha256"]),
                                "theory_manifest_sha256": np.asarray(theory["manifest"]["sha256"])})
        del kappa_alm, phi_alm
        gc.collect()
        if "nthreads" not in inspect.signature(lenspyx.alm2lenmap).parameters:
            raise RuntimeError("Installed lenspyx does not support the verified nthreads API")
        remapped = lenspyx.alm2lenmap(primary_alms, deflection_alm,
                                    geometry=("healpix", {"nside": NSIDE}),
                                    epsilon=EPSILON, verbose=1, nthreads=threads, pol=True)
        maps = np.asarray(remapped, dtype=np.float64)
        del primary_alms, deflection_alm, remapped
        gc.collect()
        if maps.shape != (3, 12 * NSIDE ** 2):
            raise ValueError(f"Unexpected lenspyx TQU shape: {maps.shape}")
        map_statistics = [full_statistics(maps[i]) for i in range(3)]
        versions = installed_versions(lenspyx)
        metadata: dict[str, Any] = {
            "schema": SCHEMA, "status": "COMPLETE", "created_utc": utc_now(),
            "realization_id": realization_id, "master_seed": master_seed, "seed": contract["seed"],
            "seed_policy": SEED_POLICY, "contract": contract,
            "contract_sha256": canonical_sha256(contract), "parameters": generation_parameters(),
            "source_kappa": source, "source_kappa_header": header,
            "source_kappa_units": "dimensionless; absent source unit header is documented rather than invented",
            "theory_archive": theory["archive"], "theory_manifest": theory["manifest"],
            "theory_configuration": theory["configuration"],
            "theory_configuration_canonical_sha256": theory["configuration_canonical_sha256"],
            "theory_normalization": theory["normalization"],
            "theory_phi_used_for_remapping": False,
            "coordinate_provenance": coordinate, "producer_source": source_records,
            "versions": versions, "finite_statistics": {"source_kappa": source_statistics,
                "generated_TQU": map_statistics},
            "input_spectrum_contract": {
                "unlensed_input_C_ell": {"shape": [4, PRIMARY_LMAX + 1], "columns": COLUMNS, "units": "uK_CMB^2 raw C_l"},
                "unlensed_sampled_C_ell": {"shape": [4, PRIMARY_LMAX + 1], "columns": COLUMNS, "units": "uK_CMB^2 raw C_l"},
                "kappa_phi_deflection": {"shape": [DEFLECTION_LMAX + 1], "units": "dimensionless raw C_l", "phi_L0": "zero", "phi_L1": "retained"},
                "phi_transfer": "kappa_to_phi: 2/[L(L+1)] for L>=1; L0=0",
                "deflection_transfer": "kappa_to_gradient_deflection: 2/sqrt[L(L+1)] for L>=1; L0=0",
                "phi_PP": "[L(L+1)]^2 C_phi_phi/(2*pi); dimensionless; L0=0, L1 retained",
                "cross_spectra": ["kappa_phi_C_ell", "kappa_deflection_C_ell", "phi_deflection_C_ell"],
                "construction": "alm2cl from existing producer alms; no additional map2alm"},
            "execution": {"host": socket.gethostname(), "pid": os.getpid(),
                "threads_requested": threads, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
                "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
                "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
                "SLURM_JOB_ID": os.environ.get("SLURM_JOB_ID"),
                "SLURM_ARRAY_TASK_ID": os.environ.get("SLURM_ARRAY_TASK_ID"),
                "mpi_rank": int(os.environ.get("OMPI_COMM_WORLD_RANK", os.environ.get("PMI_RANK", "0"))),
                "elapsed_seconds": None},
            "old_cmb_path": str(paths["old_cmb"]), "old_cmb_modified": False,
            "reused": False,
            "limitations": ["Primary and kappa cutoffs reproduce the retained producer settings, not a high-L convergence certification",
                            "The new deterministic seeds do not recover unrecorded original realizations",
                            "A uint32 seed is not a mathematical guarantee against collisions",
                            "The input kappa need not be Gaussian; a power-spectrum-only CAMB lensing reference does not reproduce all higher-order statistics"]}
        with tempfile.TemporaryDirectory(prefix=".lensedcmb_new_stage_", dir=paths["directory"]) as stage:
            staged = {k: Path(stage) / paths[k].name for k in ["fits", "spectra", "json"]}
            hp.write_map(staged["fits"], maps, dtype=np.float64, nest=False,
                         coord=coordinate["output_coordinate"], fits_IDL=True,
                         column_names=["T", "Q", "U"], column_units=["uK_CMB"] * 3,
                         extra_header=[("POLCCONV", "COSMO"), ("REALID", realization_id),
                                       ("RNGSEED", contract["seed"]), ("THSHA", theory["archive"]["sha256"]),
                                       ("KPSHA", source["sha256"]), ("UNLLMAX", PRIMARY_LMAX),
                                       ("DEFLMAX", DEFLECTION_LMAX), ("KAPITER", KAPPA_ITER),
                                       ("REMEPS", EPSILON)], overwrite=False)
            del maps
            gc.collect()
            with fits.open(staged["fits"], mode="update", memmap=True) as hdus:
                for hdu in hdus:
                    hdu.add_checksum()
            serialized_statistics = _read_output_for_verification(staged["fits"], coordinate["output_coordinate"])
            if serialized_statistics != map_statistics:
                raise RuntimeError("Float64 serialized TQU statistics differ from generated maps")
            metadata["finite_statistics"]["output_serialized_TQU"] = serialized_statistics
            np.savez_compressed(staged["spectra"], **spectrum_arrays)
            # Input identities are checked again immediately before publication.
            if file_record(paths["kappa"]) != source or file_record(theory["archive"]["path"]) != theory["archive"] or file_record(theory["manifest"]["path"]) != theory["manifest"]:
                raise RuntimeError("Scientific input changed during generation")
            backup = backup_existing(paths) if overwrite else None
            for key in ["spectra", "fits"]:
                atomic_publish(staged[key], paths[key], overwrite=overwrite)
            metadata["output"] = file_record(paths["fits"])
            metadata["input_spectra"] = file_record(paths["spectra"])
            metadata["overwrite_backup"] = backup
            metadata["execution"]["elapsed_seconds"] = time.monotonic() - started
            staged["json"].write_text(json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n")
            atomic_publish(staged["json"], paths["json"], overwrite=overwrite)
        return metadata


def parse_ids(*, ids: str | None = None, ids_file: Path | str | None = None,
              start_id: int | None = None, count: int | None = None) -> list[int]:
    modes = int(ids is not None) + int(ids_file is not None) + int(start_id is not None or count is not None)
    if modes != 1:
        raise ValueError("Specify exactly one of --ids, --ids-file, or --start-id with --count")
    if ids is not None:
        values = [int(x.strip()) for x in ids.split(",") if x.strip()]
    elif ids_file is not None:
        lines = Path(ids_file).read_text().splitlines()
        values = [int(word) for line in lines for word in line.split("#", 1)[0].replace(",", " ").split()]
    else:
        checked_integer(start_id, "start_id", maximum=99999)
        checked_integer(count, "count", minimum=1, maximum=100000)
        values = list(range(start_id, start_id + count))
    if not values or len(set(values)) != len(values):
        raise ValueError("Realization IDs must be nonempty and unique")
    for value in values:
        checked_integer(value, "realization_id", maximum=99999)
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mmdl-root", type=Path, required=True)
    parser.add_argument("--theory", type=Path, required=True)
    parser.add_argument("--theory-manifest", type=Path)
    parser.add_argument("--ids", help="Comma-separated realization IDs")
    parser.add_argument("--ids-file", type=Path, help="Whitespace/comma-separated IDs, with # comments")
    parser.add_argument("--start-id", type=int)
    parser.add_argument("--count", type=int)
    parser.add_argument("--master-seed", type=int, default=20261011)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--owner-coordinate", choices=["G", "C", "E"], default="G",
                        help="Fallback only for headerless kappa; no rotation is performed")
    parser.add_argument("--execute", action="store_true", help="Write maps; default is a read-only plan")
    parser.add_argument("--overwrite", action="store_true", help="Explicitly replace new-prefix outputs after preserving backups")
    parser.add_argument("--manifest", type=Path, help="Campaign plan/result JSON outside the source data tree")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ids = parse_ids(ids=args.ids, ids_file=args.ids_file, start_id=args.start_id, count=args.count)
    checked_integer(args.master_seed, "master_seed")
    configure_threads(args.threads)
    if args.overwrite and not args.execute:
        raise ValueError("--overwrite requires --execute")
    if args.manifest is not None:
        try:
            args.manifest.resolve().relative_to(args.mmdl_root.resolve())
        except ValueError:
            pass
        else:
            raise ValueError("Campaign manifest must be outside the mmDL source data tree")
    theory = load_theory(args.theory, args.theory_manifest)
    records = production_source_records()
    plan = {"schema": SCHEMA, "mode": "execute" if args.execute else "read-only plan",
            "created_utc": utc_now(), "ids": ids, "master_seed": args.master_seed,
            "parameters": generation_parameters(), "theory_archive": theory["archive"],
            "theory_manifest": theory["manifest"], "producer_source": records,
            "tasks": [], "errors": []}
    # No MPI dependency: independent Slurm array tasks or serial ID batches
    # call the same task seed and exact per-ID contract.
    for rid in ids:
        paths = output_paths(args.mmdl_root, rid)
        try:
            if args.execute:
                result = generate_one(args.mmdl_root, theory, rid, args.master_seed,
                                      threads=args.threads, owner_coordinate=args.owner_coordinate,
                                      overwrite=args.overwrite)
                task = {"realization_id": rid, "output": result["output"],
                        "input_spectra": result["input_spectra"], "sidecar": file_record(paths["json"]),
                        "seed": result["seed"], "reused": result["reused"]}
            else:
                header = read_kappa_header(paths["kappa"])
                coord = coordinate_provenance(header, args.owner_coordinate)
                source = file_record(paths["kappa"])
                task = {"realization_id": rid, "seed": realization_seed(args.master_seed, rid),
                        "source_kappa": source, "source_kappa_header": header,
                        "coordinate_provenance": coord,
                        "outputs": {k: str(paths[k]) for k in ["fits", "spectra", "json"]},
                        "existing_outputs": {k: paths[k].exists() for k in ["fits", "spectra", "json"]},
                        "contract_sha256": canonical_sha256(make_contract(rid, args.master_seed, source, theory, coord))}
            plan["tasks"].append(task)
            print(json.dumps(task, sort_keys=True, allow_nan=False), flush=True)
        except Exception as error:
            plan["errors"].append({"realization_id": rid, "type": type(error).__name__, "message": str(error)})
            print(json.dumps(plan["errors"][-1]), file=sys.stderr, flush=True)
            break
    plan["status"] = "FAILED" if plan["errors"] else ("COMPLETE" if args.execute else "PLAN_COMPLETE_NO_MAP_WRITES")
    if args.manifest is not None:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", prefix=".campaign_", suffix=".json", dir=args.manifest.parent, delete=False) as handle:
            json.dump(plan, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            temp = Path(handle.name)
        os.replace(temp, args.manifest)
    return int(bool(plan["errors"]))


if __name__ == "__main__":
    raise SystemExit(main())
