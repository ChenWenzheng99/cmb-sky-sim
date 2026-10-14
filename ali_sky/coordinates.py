"""Coordinate changes with the Stokes tangent-basis rotation."""
import healpy as hp
import numpy as np


def rotate_map(maps, source, target, chunk=1048576):
    """Interpolate source native pixels onto the target RING grid in chunks.

    No HEALPix downgrade. Interpolation is explicitly a map transfer, with
    Q+iU multiplied by exp(2i*angle_ref) in the target tangent basis.
    """
    if source == target:
        return maps
    if source not in ("C", "G", "E") or target not in ("C", "G", "E"):
        raise ValueError("Coordinate frames must be C, G or E")
    value = np.asarray(maps)
    scalar = value.ndim == 1
    nside = hp.npix2nside(value.shape[-1])
    rotation = hp.Rotator(coord=[source, target])
    result = np.empty_like(value, dtype=np.float64)
    for start in range(0, value.shape[-1], chunk):
        stop = min(start + chunk, value.shape[-1])
        theta_t, phi_t = hp.pix2ang(nside, np.arange(start, stop))
        theta_s, phi_s = rotation.I(theta_t, phi_t)
        if scalar:
            result[start:stop] = hp.get_interp_val(value, theta_s, phi_s)
        else:
            result[0, start:stop] = hp.get_interp_val(value[0], theta_s, phi_s)
            polarization = (hp.get_interp_val(value[1], theta_s, phi_s)
                            + 1j * hp.get_interp_val(value[2], theta_s, phi_s))
            polarization *= np.exp(2j * rotation.angle_ref(theta_s, phi_s))
            result[1, start:stop], result[2, start:stop] = polarization.real, polarization.imag
    return result


def rotate_mask(mask, source, target, chunk=1048576):
    """Nearest-pixel mask rotation preserves binary or weighted support values."""
    if source == target:
        return mask
    nside = hp.npix2nside(mask.size)
    rotation = hp.Rotator(coord=[source, target])
    result = np.empty_like(mask)
    for start in range(0, mask.size, chunk):
        stop = min(start + chunk, mask.size)
        theta, phi = hp.pix2ang(nside, np.arange(start, stop))
        theta, phi = rotation.I(theta, phi)
        result[start:stop] = mask[hp.ang2pix(nside, theta, phi)]
    return result
