"""The archived production monochromatic SEDs, in thermodynamic CMB units."""
import numpy as np


def response(component, frequency_GHz):
    frequency = float(frequency_GHz)
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("Frequency must be finite and positive in GHz")
    h, k, tcmb = 6.626068e-34, 1.3806503e-23, 2.726
    hr, kr = 6.62607004e-34, 1.38064852e-23
    if component == "ksz":
        return 1.0
    if component == "tsz":
        x = hr * frequency * 1e9 / (kr * tcmb)
        return x / np.tanh(x / 2) - 4
    if component == "cib":
        x = h * frequency * 1e9 / (k * tcmb)
        dbdt = (frequency * 1e9) ** 3 * np.exp(x) / np.expm1(x) ** 2 * x / tcmb
        return frequency ** 4.59 / np.expm1(h * frequency * 1e9 / (k * 8.35)) / dbdt
    if component == "radio":
        z = np.log(hr * frequency * 1e9 / (kr * tcmb))
        return np.exp(-2.7714227059838192 * z + .019242069853106277 * z ** 2
                      + .36813523402646997 * z ** 3)
    raise ValueError(f"Unknown extragalactic component: {component}")


def scale_factor(component, target_GHz, source_GHz=148):
    source = response(component, source_GHz)
    if component == 'tsz' and abs(source) < 1e-12:
        raise ValueError("Cannot rescale a tSZ source at its spectral null")
    if not np.isfinite(source) or source == 0:
        raise ValueError('The source SED response must be finite and nonzero')
    value = response(component, target_GHz) / source
    if not np.isfinite(value):
        raise ValueError("SED scale is not finite")
    return float(value)
