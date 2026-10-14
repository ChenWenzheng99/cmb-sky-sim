"""Reproducible Ali SAT/LAT white plus red scalar T/Q/U noise.

The numerical law and ``ell + 1`` array convention reproduce the active
``mpi_red_noise.py`` / ``noise_making.obs_noise`` generator.  This module has
no MPI dependency: a task is identified by instrument, frequency and
realization ID.  Q/U are separate scalar Gaussian draws, not spin-2 draws.
"""
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import numbers
import os
from pathlib import Path
import platform
import re
import threading

import healpy as hp
import numpy as np


_RNG_LOCK = threading.Lock()
_UINT32_MODULUS = 2 ** 32
_UINT64_MAX = 2 ** 64 - 1
_SEED_POLICY = "ali-sky-noise-seedsequence-v1"
_FREQUENCIES = (27, 40, 90, 150, 225, 280)


def _real(name, value, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a finite real number")
    value = float(value)
    if not np.isfinite(value) or (value <= 0 if positive else value < 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return value


def _integer(name, value, *, minimum=0, maximum=_UINT64_MAX):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{name} must be an integer")
    value = int(value)
    if not minimum <= value <= maximum:
        raise ValueError(f"{name} must lie in [{minimum}, {maximum}]")
    return value


def _instrument(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}", value):
        raise ValueError("instrument must contain 1-64 letters, digits, '_' or '-' and start with a letter")
    return value.upper()


def _geometry(nside, scale_factor):
    nside = _integer("nside", nside, minimum=1, maximum=2 ** 29)
    scale_factor = _integer("scale_factor", scale_factor, minimum=1, maximum=2 ** 29)
    generation_nside = nside * scale_factor
    # The legacy ud_grade call reorders through NEST, even for scale_factor=1.
    if not hp.isnsideok(nside, nest=True) or not hp.isnsideok(generation_nside, nest=True):
        raise ValueError("nside and nside*scale_factor must be valid power-of-two HEALPix NSIDEs")
    return nside, scale_factor, generation_nside


@dataclass(frozen=True)
class NoiseChannel:
    """One channel; amplitudes are in thermodynamic uK_CMB units.

    ``N_red_T`` is the temperature C_ell coefficient in uK_CMB**2.
    ``n_red_P=None`` uses the white-noise level for the polarization red term.
    The beam width is recorded for later sky simulation and is not applied to
    the generated noise. Finite positive or negative spectral slopes are
    allowed; configurations producing non-finite spectra are rejected.
    """

    frequency_GHz: float
    nlev_uK_arcmin: float
    beam_fwhm_arcmin: float
    N_red_T: float
    l_knee_T: float
    alpha_knee_T: float
    n_red_P: float | None
    l_knee_P: float
    alpha_knee_P: float
    include_red: bool = True

    def __post_init__(self):
        for name in ("frequency_GHz", "l_knee_T", "l_knee_P"):
            object.__setattr__(self, name, _real(name, getattr(self, name), positive=True))
        for name in ("nlev_uK_arcmin", "beam_fwhm_arcmin", "N_red_T"):
            object.__setattr__(self, name, _real(name, getattr(self, name)))
        if self.n_red_P is not None:
            object.__setattr__(self, "n_red_P", _real("n_red_P", self.n_red_P))
        for name in ("alpha_knee_T", "alpha_knee_P"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, numbers.Real):
                raise ValueError(f"{name} must be a finite real number")
            value = float(value)
            if not np.isfinite(value):
                raise ValueError(f"{name} must be a finite real number")
            object.__setattr__(self, name, value)
        if not isinstance(self.include_red, (bool, np.bool_)):
            raise ValueError("include_red must be a boolean")
        object.__setattr__(self, "include_red", bool(self.include_red))


def channel_from_dict(mapping):
    """Construct a channel from the exact public dataclass field names."""
    if not hasattr(mapping, "keys"):
        raise ValueError("A channel must be a mapping")
    unknown = set(mapping) - {field.name for field in fields(NoiseChannel)}
    if unknown:
        raise ValueError(f"Unknown noise channel fields: {sorted(unknown)}")
    try:
        return NoiseChannel(**dict(mapping))
    except TypeError as error:
        raise ValueError(f"Incomplete noise channel: {error}") from error


def default_channels(instrument):
    """Return the six active SAT or LAT channels in the archived generator."""
    instrument = _instrument(instrument)
    if instrument == "SAT":
        nlev = (18.27, 15.51, 3.60, 4.57, 11.46, 22.49)
        beam = (97.34, 65.70, 16.20, 9.70, 10.60, 8.70)
        red_t = (0,) * 6
        knee_p = (30, 30, 50, 50, 70, 100)
        alpha_p = (-2.4, -2.4, -2.5, -3.0, -3.0, -3.0)
    elif instrument == "LAT":
        nlev = (52, 27, 5.8, 6.3, 15, 37)
        beam = (7.4, 5.1, 2.2, 1.4, 1.0, 0.9)
        red_t = tuple(value / (8760 * 3600) for value in (100, 39, 230, 1500, 17000, 31000))
        knee_p = (700,) * 6
        alpha_p = (-1.4,) * 6
    else:
        raise ValueError("Default noise channels are available only for SAT or LAT")
    return [NoiseChannel(freq, white, fwhm, red, 1000, -3.5, white, knee, slope)
            for freq, white, fwhm, red, knee, slope in
            zip(_FREQUENCIES, nlev, beam, red_t, knee_p, alpha_p)]


def _white_noise(channel):
    try:
        white = (np.pi / (60 * 180) * channel.nlev_uK_arcmin) ** 2
    except OverflowError as error:
        raise ValueError("White noise parameters overflow the generation spectra") from error
    if not np.isfinite(white):
        raise ValueError("White noise parameters must yield a finite spectrum")
    return white


def noise_spectra(channel, nside, scale_factor=1):
    """Return legacy-indexed temperature and scalar Q/U C_ell arrays.

    Array element k evaluates the model at ell=k+1; the returned length is
    3*generation_nside-1, so synfast uses lmax=3*generation_nside-2.
    No beam or pixel window is included. At scale_factor>1 these are the
    generation spectra, before the original pixel averaging with ud_grade.
    """
    if not isinstance(channel, NoiseChannel):
        channel = channel_from_dict(channel)
    _, _, generation_nside = _geometry(nside, scale_factor)
    ell_model = np.arange(1, 3 * generation_nside)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            white = _white_noise(channel)
            if channel.include_red:
                red_p = white if channel.n_red_P is None else (np.pi / (60 * 180) * channel.n_red_P) ** 2
                nt = channel.N_red_T * (ell_model / channel.l_knee_T) ** channel.alpha_knee_T + white
                npol = 2 * (red_p * (ell_model / channel.l_knee_P) ** channel.alpha_knee_P + white)
            else:
                nt = np.full(ell_model.size, white, dtype=np.float64)
                npol = np.full(ell_model.size, 2 * white, dtype=np.float64)
    except (FloatingPointError, OverflowError) as error:
        raise ValueError("Noise parameters overflow or invalidate the generation spectra") from error
    if not np.isfinite(nt).all() or not np.isfinite(npol).all() or np.any(nt < 0) or np.any(npol < 0):
        raise ValueError("Noise spectra must be finite and nonnegative")
    return nt, npol


def _seeds(channel, rid, master_seed, instrument):
    master_seed = _integer("master_seed", master_seed)
    rid = _integer("rid", rid)
    instrument = _instrument(instrument)
    # float.hex is exact and unaffected by JSON formatting or channel order.
    identity = json.dumps([_SEED_POLICY, instrument, channel.frequency_GHz.hex()], separators=(",", ":"))
    digest = hashlib.sha256(identity.encode("ascii")).digest()
    channel_words = [int.from_bytes(digest[start:start + 4], "little") for start in range(0, 32, 4)]
    entropy = [master_seed & 0xffffffff, master_seed >> 32,
               rid & 0xffffffff, rid >> 32, *channel_words]
    seed_t = int(np.random.SeedSequence(entropy).generate_state(1, dtype=np.uint32)[0])
    return {"T": seed_t, "Q": (seed_t + 1234) % _UINT32_MODULUS,
            "U": (seed_t + 4321) % _UINT32_MODULUS}, identity, entropy


def _map_statistics(maps):
    stats = []
    for label, field in zip(("T", "Q", "U"), maps):
        finite = 0
        unseen = 0
        for start in range(0, field.size, 1048576):
            block = field[start:start + 1048576]
            finite += int(np.isfinite(block).sum())
            unseen += int(hp.mask_bad(block).sum())
        if finite != field.size or unseen:
            raise ValueError(f"Generated {label} map contains non-finite or HEALPix UNSEEN values")
        stats.append({"field": label, "pixels": int(field.size), "finite_count": finite,
                      "unseen_count": unseen, "min": float(np.min(field)), "max": float(np.max(field))})
    return stats


def _generation_contract(channel, nside, rid, master_seed, instrument, scale_factor):
    nside, scale_factor, generation_nside = _geometry(nside, scale_factor)
    seeds, identity, entropy = _seeds(channel, rid, master_seed, instrument)
    return {
        "kind": "ali_scalar_tqu_noise", "instrument": _instrument(instrument),
        "frequency_GHz": channel.frequency_GHz, "realization_id": int(rid),
        "master_seed": int(master_seed), "seeds": seeds,
        "seed_policy": _SEED_POLICY, "seed_channel_identity": identity,
        "seed_sequence_entropy_uint32": entropy,
        "seed_offsets_modulo_2_32": {"T": 0, "Q": 1234, "U": 4321},
        "seed_independent_of": ["MPI rank", "MPI size", "task order", "channel list index"],
        "channel": asdict(channel), "nside": nside, "generation_nside": generation_nside,
        "scale_factor": scale_factor, "lmax": 3 * generation_nside - 2,
        "spectrum_array_length": 3 * generation_nside - 1, "spectrum_index_convention": "array[k] evaluates model ell=k+1",
        "N_white_uK_CMB2": float(_white_noise(channel)),
        "spectra_formula": {
            "T": "N_red_T*((k+1)/l_knee_T)**alpha_knee_T + N_white",
            "Q_U": "2*(N_red_P*((k+1)/l_knee_P)**alpha_knee_P + N_white)",
            "N_red_P": "(pi*n_red_P/10800)**2; None uses N_white; red terms zero if include_red=False"},
        "scalar_fields": ["T", "Q", "U"], "pol_argument": False,
        "spin2_polarization_generation": False, "beam_applied": False,
        "mask_applied": False, "pixel_window_applied": False, "monopole_removed": False,
        "resampling": "legacy healpy.ud_grade RING pixel averaging" if scale_factor > 1 else "legacy healpy.ud_grade at unchanged NSIDE",
        "unit": "uK_CMB", "ordering": "RING", "dtype": "float64",
        "source_convention": {"driver": "mpi_red_noise.py active SAT/LAT sections", "generator": "noise_making.py:obs_noise"},
    }


def generate_noise(channel, nside, rid, master_seed, instrument="custom", scale_factor=1):
    """Generate one float64 RING T/Q/U realization and JSON-safe provenance.

    The three legacy np.random seeds are used under a process-local lock and
    the caller's global RNG state is restored, also when generation fails.
    Channel list position and MPI rank never enter the seed derivation.
    """
    if not isinstance(channel, NoiseChannel):
        channel = channel_from_dict(channel)
    metadata = _generation_contract(channel, nside, rid, master_seed, instrument, scale_factor)
    nt, npol = noise_spectra(channel, nside, scale_factor)
    maps = np.empty((3, hp.nside2npix(metadata["nside"])), dtype=np.float64)
    with _RNG_LOCK:
        caller_state = np.random.get_state()
        try:
            for index, (label, spectrum) in enumerate(zip(("T", "Q", "U"), (nt, npol, npol))):
                np.random.seed(metadata["seeds"][label])
                realization = hp.synfast(spectrum, metadata["generation_nside"], pol=False, new=True)
                maps[index] = hp.ud_grade(realization, metadata["nside"])
        finally:
            np.random.set_state(caller_state)
    metadata.update({
        "statistics": _map_statistics(maps),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "healpy": hp.__version__},
        "thread_environment": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
    })
    return maps, metadata


def _resume(path, expected, coord):
    """Accept an existing product only if its contract and bytes are verified."""
    from .io import read_fits

    try:
        metadata = json.loads(path.with_suffix(".json").read_text())
        if not isinstance(metadata, dict):
            raise ValueError("sidecar must be a JSON object")
        for key, value in expected.items():
            if metadata.get(key) != value:
                raise ValueError(f"sidecar mismatch in {key}")
        if metadata.get("coordinate") != coord or metadata.get("polarization_convention") != "COSMO":
            raise ValueError("sidecar coordinate or polarization convention mismatch")
        output = metadata["output"]
        if not isinstance(output, dict):
            raise ValueError("sidecar output record must be an object")
        maps, actual_coord, info = read_fits(path, expected["nside"])
        if actual_coord != coord or output.get("path") != str(path.resolve()):
            raise ValueError("FITS coordinate or output path mismatch")
        if output.get("sha256") != info["sha256"] or output.get("bytes") != info["bytes"]:
            raise ValueError("output hash or size mismatch")
        if metadata.get("statistics") != _map_statistics(maps):
            raise ValueError("map statistics mismatch")
    except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
        raise FileExistsError(f"Existing output cannot be safely reused: {path}: {error}") from error
    return {**metadata, "reused": True, "resume_validation": "contract, seeds, FITS metadata, full finite statistics and SHA256 matched"}


def generate_noisefile(config, task):
    """Generate one normalized-config task and write FITS plus JSON via IO.

    ``task=(channel_index, realization_id)``; the output layout is
    output_dir/instrument/frequency/noise_{id:05d}.fits plus a JSON sidecar.
    The configured coordinate names the synthesis grid; no rotation is made.
    """
    from .config import frequency_tag
    from .io import write_tqu

    if not isinstance(task, (tuple, list)) or len(task) != 2:
        raise ValueError("A noise task must be (channel_index, realization_id)")
    index = _integer("channel_index", task[0])
    rid = _integer("rid", task[1])
    channels = config["channels"]
    if index >= len(channels):
        raise ValueError("channel_index is outside the configured channels")
    ids = [_integer("ids entry", value) for value in config["ids"]]
    if rid not in ids:
        raise ValueError("realization_id is outside the configured ids")
    channel = channels[index]
    if not isinstance(channel, NoiseChannel):
        channel = channel_from_dict(channel)
    instrument = _instrument(config["instrument"])
    coord = config["coordinate"]
    if coord not in ("C", "G", "E"):
        raise ValueError("coordinate must be C, G or E")
    overwrite = config["overwrite"]
    if not isinstance(overwrite, bool):
        raise ValueError("overwrite must be a boolean")
    token = frequency_tag(channel.frequency_GHz)
    path = Path(config["output_dir"]) / instrument / token / f"noise_{rid:05d}.fits"
    expected = _generation_contract(channel, config["nside"], rid, config["seed"], instrument, config["scale_factor"])
    if not overwrite and (path.exists() or path.with_suffix(".json").exists()):
        return _resume(path, expected, coord)
    maps, metadata = generate_noise(channel, config["nside"], rid, config["seed"],
                                    instrument, config["scale_factor"])
    metadata["coordinate"] = coord
    metadata["coordinate_assignment"] = "independent scalar fields generated directly on this coordinate-labelled RING grid; no map rotation"
    return {**write_tqu(path, maps, coord, metadata, overwrite=overwrite), "reused": False}
