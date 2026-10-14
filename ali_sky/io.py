"""Strict native HEALPix input/output; no private pipeline dependencies."""
import hashlib
import json
import os
from pathlib import Path
import tempfile

import healpy as hp
import numpy as np
from astropy.io import fits

from . import __version__


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_info(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "sha256": sha256(path)}


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def check_map(maps, nside, fields=3):
    maps = np.asarray(maps)
    expected = (3, hp.nside2npix(nside)) if fields == 3 else (hp.nside2npix(nside),)
    if maps.shape != expected:
        raise ValueError(f"Map shape {maps.shape}; expected {expected}")
    for field in np.atleast_2d(maps):
        for start in range(0, field.size, 1048576):
            block = field[start:start + 1048576]
            if not np.isfinite(block).all() or hp.mask_bad(block).any():
                raise ValueError("Map contains non-finite or HEALPix UNSEEN pixels")
    return maps


def fits_frame(path, assume_coord=None, require_unit=True, assume_unit=None):
    header = fits.getheader(path, 1)
    coord = str(header.get("COORDSYS", assume_coord or "")).strip().upper()
    if coord not in ("C", "G", "E"):
        raise ValueError(f"Missing/invalid COORDSYS in {path}; supply an explicit source frame")
    units = [str(header.get(f"TUNIT{i}", assume_unit or "")).strip()
             for i in range(1, header.get("TFIELDS", 1) + 1)]
    if require_unit and any(unit not in ("uK_CMB", "microK_CMB") for unit in units):
        raise ValueError(f"Expected thermodynamic uK_CMB units in {path}, found {units}")
    pol = str(header.get("POLCCONV", "COSMO")).strip().upper()
    if pol != "COSMO":
        raise ValueError(f"Unsupported polarization convention {pol} in {path}")
    return coord, dict(header)


def read_fits(path, nside, fields=3, assume_coord=None, require_unit=True, assume_unit=None):
    coord, header = fits_frame(path, assume_coord, require_unit, assume_unit)
    if header.get("NSIDE") != nside:
        raise ValueError(f"{path}: NSIDE={header.get('NSIDE')}, requested {nside}; no automatic resampling")
    maps = hp.read_map(str(path), field=(0, 1, 2) if fields == 3 else 0,
                       dtype=np.float64, nest=False)
    check_map(maps, nside, fields)
    info = file_info(path)
    info.update({'coordinate': coord, 'coordinate_from_FITS_header': 'COORDSYS' in header,
                 'assumed_coordinate_if_header_missing': assume_coord,
                 'unit_assumption_if_header_missing': assume_unit})
    return maps, coord, info


def read_noise(path, nside, assume_coord=None):
    path = Path(path)
    if path.suffix.lower() in (".fits", ".fit", ".gz"):
        return read_fits(path, nside, assume_coord=assume_coord, assume_unit=None)
    if path.suffix.lower() != ".npy":
        raise ValueError("Noise input must be T/Q/U FITS or NPY")
    sidecar = path.with_suffix(".json")
    meta = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    coord = meta.get("coordinate", assume_coord)
    if coord not in ("C", "G", "E"):
        raise ValueError("A coordinate frame is required for NPY noise")
    data = np.load(path, mmap_mode="r", allow_pickle=False)
    npix = hp.nside2npix(nside)
    if data.shape == (3, npix):
        maps = np.array(data, dtype=np.float64)
    elif data.shape == (3, 2, npix):
        maps = np.empty((3, npix), dtype=np.float64)
        for stokes in range(3):
            for start in range(0, npix, 1048576):
                stop = min(start + 1048576, npix)
                if not np.array_equal(data[stokes, 0, start:stop], np.arange(start, stop)):
                    raise ValueError("Legacy noise indices must cover the complete monotone RING grid")
            maps[stokes] = data[stokes, 1]
    else:
        raise ValueError(f"Invalid noise shape {data.shape}; expected (3,npix) or (3,2,npix)")
    check_map(maps, nside)
    return maps, coord, file_info(path)


def write_tqu(path, maps, coord, metadata, overwrite=False):
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists: {path}; enable overwrite explicitly")
    if coord not in ("C", "G", "E"):
        raise ValueError("coordinate must be C, G or E")
    nside = hp.npix2nside(np.asarray(maps).shape[-1])
    check_map(maps, nside)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".fits", dir=path.parent)
    os.close(descriptor)
    try:
        extra = [("POLCCONV", "COSMO"), ("CREATOR", f"ali_sky {__version__}")]
        if "frequency_GHz" in metadata:
            extra.append(("FREQGHZ", float(metadata["frequency_GHz"])))
        hp.write_map(temporary, maps, dtype=np.float64, coord=coord, nest=False,
                     column_names=["T", "Q", "U"], column_units=["uK_CMB"] * 3,
                     extra_header=extra, overwrite=True)
        with fits.open(temporary, mode="update", memmap=True) as hdul:
            for hdu in hdul:
                hdu.add_checksum()
        if overwrite:
            os.replace(temporary, path)
        else:
            # Same-directory hard-link creation is atomic and cannot clobber a
            # concurrent job's output after the initial existence check.
            os.link(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    meta = {**metadata, "coordinate": coord, "nside": nside, "ordering": "RING",
            "polarization_convention": "COSMO", "unit": "uK_CMB", "output": file_info(path)}
    atomic_json(path.with_suffix(".json"), meta)
    return meta
