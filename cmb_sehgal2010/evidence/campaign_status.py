#!/usr/bin/env python3
"""Read-only source inventory and bounded-I/O CMB campaign status checks.

capture uses stat only: it never hashes old CMB or kappa pixel data.
status verifies every old size/mtime_ns/inode tuple, COMPLETE new sidecars,
input/source/seed/producer contracts, new file stat, lightweight FITS headers,
and hashes of JSON and small input_spectra archives. It does not independently
rehash new full-resolution maps or recompute FITS data checksums. Full serialized
map checks are the frozen producer's recorded responsibility; independent
native pilot validation remains a separate evidence layer.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import socket
import struct
import sys
import tempfile


SCHEMA = "sehgal2010_cmb_campaign_status_v1"
NPIX = 12 * 2048 ** 2
STAT_FIELDS = ["bytes", "mtime_ns", "inode"]
SEED_DOMAIN = "sehgal2010-primary-mmdl-kappa-native2048-v1"
SEED_POLICY = "numpy.SeedSequence(domain_sha256_u32, master_seed_u64, realization_id_u64); uint32_v1"
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def stat_record(path):
    path = Path(path)
    info = path.stat()
    if not path.is_file() or path.is_symlink():
        raise ValueError(f"Expected a regular non-symlink file: {path}")
    return {"path": str(path.resolve()), "bytes": info.st_size,
            "mtime_ns": info.st_mtime_ns, "inode": info.st_ino,
            "device": info.st_dev, "mode": info.st_mode & 0o777}


def small_file_record(path):
    before = stat_record(path)
    record = {**before, "sha256": sha256(path)}
    after = stat_record(path)
    if before != after:
        raise RuntimeError(f"Small artifact changed while hashing: {path}")
    return record


def tool_record(external_sha256=None):
    path = Path(__file__)
    if path.is_file():
        return small_file_record(path)
    if external_sha256 and SHA_RE.fullmatch(external_sha256):
        return {"path": "stdin", "sha256": external_sha256,
                "identity": "Hash of the exact local helper supplied by the SSH caller"}
    return {"path": str(path), "sha256": None}


def paths_for(root, rid):
    directory = Path(root) / f"{rid:05d}"
    prefix = f"lensedcmb_new_{rid:05d}"
    return {"directory": directory, "old": directory / f"lensedcmb_{rid:05d}.fits",
            "kappa": directory / f"kappa_{rid:05d}.fits", "fits": directory / (prefix + ".fits"),
            "json": directory / (prefix + ".json"), "spectra": directory / (prefix + ".input_spectra.npz"),
            "lock": directory / (prefix + ".lock")}


def numeric_directory_ids(root):
    return sorted(int(entry.name) for entry in Path(root).iterdir()
                  if entry.is_dir() and re.fullmatch(r"[0-9]{5}", entry.name))


def capture_inventory(root, ids=None, external_tool_sha256=None):
    root = Path(root).resolve(strict=True)
    ids = list(range(500)) if ids is None else list(ids)
    observed = numeric_directory_ids(root)
    if observed != sorted(ids):
        raise ValueError("Numeric source directories do not exactly match the expected realization IDs")
    rows = []
    old_count = kappa_count = new_count = 0
    unexpected_source_filenames = []
    for rid in ids:
        paths = paths_for(root, rid)
        old, kappa = stat_record(paths["old"]), stat_record(paths["kappa"])
        for entry in paths["directory"].iterdir():
            if entry.is_file() and re.fullmatch(r"lensedcmb_[0-9]{5}\.fits", entry.name):
                old_count += 1
                if entry != paths["old"]:
                    unexpected_source_filenames.append(str(entry))
            if entry.is_file() and re.fullmatch(r"kappa_[0-9]{5}\.fits", entry.name):
                kappa_count += 1
                if entry != paths["kappa"]:
                    unexpected_source_filenames.append(str(entry))
        new_count += int(paths["fits"].exists())
        rows.append({"realization_id": rid, "old_cmb": old, "kappa": kappa,
                     "new_outputs_present": {key: paths[key].exists() for key in ["fits", "json", "spectra", "lock"]}})
    return {"schema": SCHEMA, "kind": "BEFORE_SOURCE_STAT_INVENTORY", "created_utc": utc_now(),
            "host": socket.gethostname(), "root": str(root), "ids": ids,
            "directory_ids": observed, "old_cmb_count": old_count, "kappa_count": kappa_count,
            "source_filename_counts_exact": old_count == len(ids) and kappa_count == len(ids) and not unexpected_source_filenames,
            "unexpected_source_filenames": unexpected_source_filenames,
            "existing_new_fits_count": new_count, "realizations": rows,
            "tool": tool_record(external_tool_sha256), "python_version": platform.python_version(),
            "large_source_files_hashed": False,
            "comparison_fields": STAT_FIELDS,
            "scope": "Baseline size, mtime_ns and inode identity; this is not a full-content SHA256 audit"}


def independent_seed(master_seed, rid):
    import numpy as np
    if isinstance(master_seed, bool) or not isinstance(master_seed, int) or not 0 <= master_seed < 2 ** 64:
        raise ValueError("master_seed must be an unsigned 64-bit integer")
    words = list(struct.unpack("<8I", hashlib.sha256(SEED_DOMAIN.encode()).digest()))
    words += [master_seed & 0xffffffff, master_seed >> 32, rid & 0xffffffff, rid >> 32]
    return int(np.random.SeedSequence(words).generate_state(1, dtype=np.uint32)[0])


def load_producer(path):
    path = Path(path).resolve(strict=True)
    spec = importlib.util.spec_from_file_location("campaign_expected_producer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, small_file_record(path)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def header_record(path):
    """Read header cards only; never access hdu.data or verify data checksums."""
    from astropy.io import fits
    keys = ["NSIDE", "ORDERING", "COORDSYS", "POLCCONV", "TFIELDS", "INDXSCHM", "NAXIS1", "NAXIS2",
            "FIRSTPIX", "LASTPIX", "REALID", "RNGSEED", "THSHA", "KPSHA",
            "UNLLMAX", "DEFLMAX", "KAPITER", "REMEPS", "TTYPE1", "TTYPE2", "TTYPE3",
            "TUNIT1", "TUNIT2", "TUNIT3", "TFORM1", "TFORM2", "TFORM3"]
    with fits.open(path, memmap=True, checksum=False, lazy_load_hdus=True) as hdus:
        result = {key: hdus[1].header.get(key) for key in keys}
        result["checksum_cards"] = [{"CHECKSUM": hdu.header.get("CHECKSUM"),
                                     "DATASUM": hdu.header.get("DATASUM")} for hdu in hdus]
    result["pixel_data_read"] = False
    result["data_checksums_recomputed"] = False
    return result


def verify_stat_identity(current, reference, description):
    for key in STAT_FIELDS:
        require(current.get(key) == reference.get(key), f"{description}: {key} changed")


def verify_finite_statistics(metadata):
    saved = metadata.get("finite_statistics", {})
    source = saved.get("source_kappa", {})
    require(source.get("shape") == [NPIX] and source.get("count") == NPIX and
            source.get("finite_count") == NPIX and source.get("unseen_count") == 0,
            "Full-array source finite/UNSEEN record is incomplete")
    generated = saved.get("generated_TQU")
    serialized = saved.get("output_serialized_TQU")
    require(isinstance(generated, list) and len(generated) == 3 and generated == serialized,
            "Serialized full-array finite records differ from generated TQU")
    for name, field in [("source_kappa", source)] + [(f"serialized_{i}", f) for i, f in enumerate(serialized)]:
        require(field.get("shape") == [NPIX] and field.get("count") == NPIX and
                field.get("finite_count") == NPIX and field.get("unseen_count") == 0 and
                field.get("dtype") == "float64", f"Invalid full-array finite record: {name}")
        require(all(isinstance(field.get(k), (float, int)) and math.isfinite(field[k])
                    for k in ["min", "max", "mean", "population_std"]), f"Non-finite recorded summary: {name}")
        require(field["min"] <= field["max"] and field["population_std"] >= 0,
                f"Invalid recorded extrema/standard deviation: {name}")


def verify_header(header, metadata, rid, expected_seed, producer):
    expected = {"NSIDE": 2048, "ORDERING": "RING", "COORDSYS": metadata["coordinate_provenance"]["output_coordinate"],
                "POLCCONV": "COSMO", "TFIELDS": 3, "INDXSCHM": "IMPLICIT", "FIRSTPIX": 0,
                "LASTPIX": NPIX - 1, "REALID": rid, "RNGSEED": expected_seed,
                "THSHA": metadata["theory_archive"]["sha256"], "KPSHA": metadata["source_kappa"]["sha256"],
                "UNLLMAX": producer.PRIMARY_LMAX, "DEFLMAX": producer.DEFLECTION_LMAX,
                "KAPITER": producer.KAPPA_ITER, "REMEPS": producer.EPSILON}
    require(all(header.get(k) == v for k, v in expected.items()), "New FITS header differs from producer contract")
    require([header.get(f"TTYPE{i}") for i in [1, 2, 3]] == ["T", "Q", "U"], "New FITS column order is not TQU")
    require(all(header.get(f"TUNIT{i}") == "uK_CMB" and str(header.get(f"TFORM{i}", "")).endswith("D")
                for i in [1, 2, 3]), "New FITS dtype/units are not float64 uK_CMB")
    repeats = [re.fullmatch(r"([0-9]*)D", str(header.get(f"TFORM{i}", ""))) for i in [1, 2, 3]]
    require(all(match is not None for match in repeats), "Unexpected floating-point table format")
    repeat_values = [int(match.group(1) or "1") for match in repeats]
    require(all(repeat * header.get("NAXIS2", 0) == NPIX for repeat in repeat_values) and
            header.get("NAXIS1") == 8 * sum(repeat_values), "FITS table dimensions do not contain native full-sky TQU")
    require(all(isinstance(cards.get("CHECKSUM"), str) and len(cards["CHECKSUM"]) == 16 and
                cards.get("DATASUM") is not None for cards in header["checksum_cards"]),
            "New FITS checksum cards are absent")


def verify_complete(paths, baseline, producer, producer_record, theory_record,
                    theory_manifest_record, master_seed, owner_coordinate):
    rid = baseline["realization_id"]
    metadata = json.loads(paths["json"].read_text())
    require(metadata.get("schema") == producer.SCHEMA and metadata.get("status") == "COMPLETE",
            "Sidecar is not a COMPLETE producer record")
    require(metadata.get("realization_id") == rid and metadata.get("master_seed") == master_seed,
            "Sidecar realization/master seed differs from the campaign")
    seed = independent_seed(master_seed, rid)
    require(metadata.get("seed") == seed and metadata.get("seed_policy") == SEED_POLICY,
            "Sidecar seed does not follow the independently derived task policy")
    require(metadata.get("parameters") == producer.generation_parameters(), "Producer parameters changed")
    require(metadata.get("old_cmb_modified") is False and metadata.get("old_cmb_path") == str(paths["old"]),
            "Old-map preservation provenance is missing or points at another file")
    require(metadata.get("theory_phi_used_for_remapping") is False, "CAMB phi was substituted for actual kappa")
    source_header = producer.read_kappa_header(paths["kappa"])
    coord = producer.coordinate_provenance(source_header, owner_coordinate)
    require(metadata.get("source_kappa_header") == source_header and metadata.get("coordinate_provenance") == coord,
            "Source header/frame provenance differs from the actual source header")
    source = metadata.get("source_kappa", {})
    require(source.get("path") == str(paths["kappa"]) and SHA_RE.fullmatch(str(source.get("sha256", ""))),
            "Source kappa identity/hash is missing")
    for key in ["bytes", "mtime_ns"]:
        require(source.get(key) == baseline["kappa"].get(key), f"Source kappa recorded {key} differs from baseline")
    for name, expected in [("theory_archive", theory_record), ("theory_manifest", theory_manifest_record)]:
        saved = metadata.get(name, {})
        require(saved.get("sha256") == expected["sha256"] and saved.get("bytes") == expected["bytes"],
                f"Recorded {name} content identity differs from the supplied theory")
    contract = producer.make_contract(rid, master_seed, source,
                                     {"archive": theory_record, "manifest": theory_manifest_record}, coord)
    require(metadata.get("contract") == contract and metadata.get("contract_sha256") == canonical_sha256(contract),
            "Complete sidecar contract/hash mismatch")
    source_code = metadata.get("producer_source", {})
    require(source_code.get("release_producer", {}).get("sha256") == producer_record["sha256"] and
            source_code.get("original_snapshot", {}).get("sha256") == producer.ORIGINAL_SHA256,
            "Recorded producer/original code hash differs from the frozen release")
    require(metadata.get("overwrite_backup") is None, "Initial campaign unexpectedly used overwrite")
    verify_finite_statistics(metadata)
    output_stat = stat_record(paths["fits"])
    saved_output = metadata.get("output", {})
    require(saved_output.get("path") == str(paths["fits"]) and SHA_RE.fullmatch(str(saved_output.get("sha256", ""))),
            "Recorded output path/hash is missing")
    for key in ["bytes", "mtime_ns"]:
        require(saved_output.get(key) == output_stat[key], f"New FITS {key} differs from produced stat")
    spectra = small_file_record(paths["spectra"])
    saved_spectra = metadata.get("input_spectra", {})
    require(saved_spectra.get("path") == str(paths["spectra"]) and all(saved_spectra.get(k) == spectra[k]
                for k in ["bytes", "mtime_ns", "sha256"]), "Small input_spectra content/stat differs from the sidecar")
    header = header_record(paths["fits"])
    verify_header(header, metadata, rid, seed, producer)
    return {"realization_id": rid, "state": "COMPLETE_STAT_HEADER_VERIFIED", "seed": seed,
            "sidecar": small_file_record(paths["json"]), "output_current_stat": output_stat,
            "output_sha256_recorded_by_producer": saved_output["sha256"],
            "output_full_sha256_recomputed_by_this_tool": False,
            "input_spectra": spectra, "source_kappa_sha256_recorded_by_producer": source["sha256"],
            "source_kappa_sha256_recomputed_by_this_tool": False,
            "producer_sha256": producer_record["sha256"], "theory_archive_sha256": theory_record["sha256"],
            "coordinate_provenance": coord, "header": header, "versions": metadata.get("versions"),
            "execution": metadata.get("execution"),
            "recorded_full_array_count_TQU": 3 * NPIX,
            "producer_full_array_serialization_checks_recorded": True,
            "stat_header_check_is_full_content_validation": False}


def campaign_status(inventory_path, producer_path, theory_path, theory_manifest_path=None,
                    master_seed=20261011, owner_coordinate="G"):
    inventory_path = Path(inventory_path).resolve(strict=True)
    baseline = json.loads(inventory_path.read_text())
    require(baseline.get("schema") == SCHEMA and baseline.get("kind") == "BEFORE_SOURCE_STAT_INVENTORY",
            "Invalid before-inventory schema")
    root = Path(baseline["root"]).resolve(strict=True)
    ids = baseline["ids"]
    producer, producer_record = load_producer(producer_path)
    # The frozen loader verifies raw-Cl conventions, array shapes and theory hashes.
    theory = producer.load_theory(theory_path, theory_manifest_path)
    theory_record, theory_manifest_record = theory["archive"], theory["manifest"]
    require(len(ids) == 500 and ids == list(range(500)), "Campaign status requires the complete 0..499 baseline")
    require(baseline.get("source_filename_counts_exact") is True and baseline.get("old_cmb_count") == 500 and
            baseline.get("kappa_count") == 500, "Before-inventory source filename counts were not exact")
    expected_seeds = [independent_seed(master_seed, rid) for rid in ids]
    require(len(set(expected_seeds)) == len(ids), "Expected uint32 seeds collide; choose a different master_seed before execution")
    require([row.get("realization_id") for row in baseline["realizations"]] == ids,
            "Before-inventory realization rows are incomplete or reordered")
    report = {"schema": SCHEMA, "kind": "CAMPAIGN_STATUS", "created_utc": utc_now(),
              "host": socket.gethostname(), "root": str(root), "ids": ids,
              "before_inventory": small_file_record(inventory_path), "producer": producer_record,
              "theory_archive": theory_record, "theory_manifest": theory_manifest_record,
              "master_seed": master_seed, "tool": tool_record(), "rows": [], "errors": [],
              "old_stat_unchanged_count": 0, "kappa_stat_unchanged_count": 0,
              "expected_seed_unique_count": len(set(expected_seeds))}
    if numeric_directory_ids(root) != baseline["directory_ids"]:
        report["errors"].append({"scope": "source directory IDs", "message": "Numeric realization directory IDs changed"})
    for row in baseline["realizations"]:
        rid = row["realization_id"]
        paths = paths_for(root, rid)
        source_errors = []
        for key, name, countkey in [("old", "old_cmb", "old_stat_unchanged_count"),
                                     ("kappa", "kappa", "kappa_stat_unchanged_count")]:
            try:
                verify_stat_identity(stat_record(paths[key]), row[name], name)
                report[countkey] += 1
            except Exception as error:
                source_errors.append({"scope": name, "type": type(error).__name__, "message": str(error)})
        present = {key: paths[key].exists() for key in ["fits", "json", "spectra", "lock"]}
        try:
            if source_errors:
                raise ValueError("Source stat identity changed: " + json.dumps(source_errors))
            if all(present[key] for key in ["fits", "json", "spectra"]):
                result = verify_complete(paths, row, producer, producer_record, theory_record,
                                         theory_manifest_record, master_seed, owner_coordinate)
                result["lock_present"] = present["lock"]
                if present["lock"]:
                    result["state"] = "COMPLETE_WITH_ACTIVE_OR_STALE_LOCK"
            else:
                result = {"realization_id": rid, "state": "IN_PROGRESS" if present["lock"] else
                          ("PARTIAL_NO_COMMIT" if any(present.values()) else "NOT_STARTED"),
                          "outputs_present": present}
        except Exception as error:
            result = {"realization_id": rid, "state": "INVALID", "error_type": type(error).__name__,
                      "message": str(error), "outputs_present": present}
            report["errors"].append({"realization_id": rid, **result})
        report["rows"].append(result)
    states = sorted({row["state"] for row in report["rows"]})
    report["state_counts"] = {state: sum(row["state"] == state for row in report["rows"]) for state in states}
    completed = [row for row in report["rows"] if row["state"] == "COMPLETE_STAT_HEADER_VERIFIED"]
    report["completed_count"] = len(completed)
    report["recorded_completed_seed_unique_count"] = len({row["seed"] for row in completed})
    require(report["recorded_completed_seed_unique_count"] == len(completed), "Completed-sidecar seeds are not unique")
    report["all_500_old_stat_tuples_unchanged"] = report["old_stat_unchanged_count"] == 500
    report["all_500_kappa_stat_tuples_unchanged"] = report["kappa_stat_unchanged_count"] == 500
    report["campaign_complete"] = len(completed) == 500 and not report["errors"]
    report["status"] = "COMPLETE_STAT_HEADER_AUDIT" if report["campaign_complete"] else (
        "INVALID" if report["errors"] else "INCOMPLETE")
    report["scope"] = {"old_maps": "All 500 size/mtime_ns/inode tuples compared; no old full-content hashes",
                       "source_kappa": "Stat identity and producer-recorded SHA256 checked; no independent full-kappa rehash",
                       "new_maps": "Stat identity and header/checksum-card presence checked; no new full-map hash or data-checksum recomputation",
                       "small_artifacts": "Every COMPLETE sidecar and input_spectra NPZ independently SHA256 hashed",
                       "pixel_validation": "Full finite/statistics/checksum validation is recorded by the frozen producer; independent pilot science evidence is separate",
                       "science_certification": False}
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    capture = sub.add_parser("capture", help="Stat-only before inventory of IDs 0..499")
    capture.add_argument("--root", type=Path, required=True)
    capture.add_argument("--external-tool-sha256", help="Exact supplied helper SHA when executed from SSH stdin")
    status = sub.add_parser("status", help="All-ID aggregate status after or during the campaign")
    status.add_argument("--inventory", type=Path, required=True)
    status.add_argument("--producer", type=Path, required=True)
    status.add_argument("--theory", type=Path, required=True)
    status.add_argument("--theory-manifest", type=Path)
    status.add_argument("--master-seed", type=int, default=20261011)
    status.add_argument("--owner-coordinate", choices=["G", "C", "E"], default="G")
    status.add_argument("--require-complete", action="store_true", help="Return nonzero unless all 500 strict checks pass")
    for item in [capture, status]:
        item.add_argument("--output", type=Path, help="Evidence JSON outside source data tree; default stdout")
    args = parser.parse_args(argv)
    if args.command == "capture":
        report = capture_inventory(args.root, external_tool_sha256=args.external_tool_sha256)
    else:
        report = campaign_status(args.inventory, args.producer, args.theory, args.theory_manifest,
                                 args.master_seed, args.owner_coordinate)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        sys.stdout.write(payload)
    else:
        try:
            args.output.resolve().relative_to(Path(report["root"]))
        except ValueError:
            pass
        else:
            raise ValueError("Status evidence output must be outside the source data tree")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile("w", dir=args.output.parent, prefix=".status_", delete=False) as handle:
            handle.write(payload)
            staged = Path(handle.name)
        os.replace(staged, args.output)
        print(json.dumps({"output": str(args.output), "sha256": sha256(args.output),
                          "status": report.get("status", report.get("kind")),
                          "completed_count": report.get("completed_count")}, sort_keys=True))
    return int(args.command == "status" and (bool(report["errors"]) or
               (args.require_complete and not report["campaign_complete"])))


if __name__ == "__main__":
    raise SystemExit(main())
