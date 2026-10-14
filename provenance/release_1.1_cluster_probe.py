"""Read-only source/runtime checks and published CMB header provenance.

This probe uses fresh tiny unit-test fixtures; it does not remap native sky
inputs or recompute their full-file hashes.
"""
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

from astropy.io import fits

from ali_sky import __version__, config


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    assert __version__ == "1.1.0"
    command = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"]
    result = subprocess.run(command, cwd=root, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    log = root / "provenance/release_1.1_cluster_tests.log"
    log.write_text(result.stdout)
    if result.returncode or "Ran 32 tests" not in result.stdout or not result.stdout.rstrip().endswith("OK"):
        raise RuntimeError(result.stdout)
    examples = []
    for path in sorted((root / "examples").glob("*.json")):
        mode = "noise" if path.name.startswith("noise_") else "observe"
        normalized = config.load_config(path, mode)
        examples.append({"path": path.relative_to(root).as_posix(), "mode": mode,
                         "frequencies_GHz": [c["frequency_GHz"] for c in normalized["channels"]],
                         "output_dir": normalized["output_dir"], "coordinate": normalized["coordinate"]})
    assert len(examples) == 6
    cfg = config.load_config(root / "examples/observe_lat150.json", "observe")
    folder = Path(cfg["repository"]["sky_dir"]) / "00000"
    current, legacy = folder / "lensedcmb_new_00000.fits", folder / "lensedcmb_00000.fits"
    assert current.is_file() and legacy.is_file()
    header = fits.getheader(current, 1)
    values = {key: header[key] for key in ("NSIDE", "ORDERING", "COORDSYS", "POLCCONV",
              "TUNIT1", "TUNIT2", "TUNIT3", "REALID", "RNGSEED", "UNLLMAX", "DEFLMAX", "THSHA", "KPSHA")}
    assert values["NSIDE"] == 2048 and values["COORDSYS"] == "G"
    assert values["UNLLMAX"] == 5999 and values["DEFLMAX"] == 6143
    assert all(values[f"TUNIT{i}"] == "uK_CMB" for i in (1, 2, 3))
    sidecar = current.with_suffix(".json")
    producer = json.loads(sidecar.read_text())
    assert producer["status"] == "COMPLETE"
    assert producer["seed"] == values["RNGSEED"]
    assert producer["theory_archive"]["sha256"] == values["THSHA"]
    assert producer["source_kappa"]["sha256"] == values["KPSHA"]
    from mpi4py import MPI
    report = {
        "status": "PASS_SCOPED_CURRENT_SOURCE_RUNTIME_CHECKS", "release": "ali_sky 1.1.0",
        "source_root": str(root), "python": sys.version, "python_executable": sys.executable,
        "runtime_dependency_reuse": "Existing immutable 1.0 runtime, read-only; no new environment relocation claimed",
        "host": platform.node(), "versions": {n: metadata.version(n) for n in
            ("numpy", "healpy", "astropy", "scipy", "pysm3", "mpi4py")},
        "MPI_vendor": MPI.get_vendor(), "test_count": 32, "test_command": command,
        "test_exit_code": result.returncode, "test_log": {"path": log.relative_to(root).as_posix(), "sha256": digest(log)},
        "example_configurations": examples,
        "tested_sources": [{"path": p.relative_to(root).as_posix(), "sha256": digest(p)}
            for p in sorted((root / "ali_sky").glob("*.py")) + sorted((root / "tests").glob("*.py"))],
        "default_CMB_provenance": {"selected_path": str(current), "legacy_path_present_but_ignored": str(legacy),
            "header": values, "sidecar": {"path": str(sidecar), "sha256": digest(sidecar)},
            "pixel_data_read": False, "full_FITS_sha256_recomputed": False,
            "primary": "Sehgal2010 CAMB primary, r=0", "lensing": "Corresponding actual mmDL kappa",
            "frame_caveat": producer["coordinate_provenance"]["caveat"]},
        "scope_limits": ["Complete software suite and example normalization; no native observation production rerun.",
            "Published CMB header and sidecar checks do not independently certify the headerless source kappa physical frame."]}
    out = root / "provenance/release_1.1_cluster_validation.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "test_count": 32, "report": str(out)}, indent=2))


if __name__ == "__main__":
    main()
