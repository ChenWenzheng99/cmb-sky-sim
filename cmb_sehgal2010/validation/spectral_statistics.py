"""Predeclared pilot statistics with stable signed-TE residuals.

Sample standard deviation describes realization scatter. SEM is used only
for the named pilot-mean consistency diagnostics. The BB diagnostics do not
pretend Gaussian diagonal theory covariance is its full sampling covariance.
"""
import numpy as np
from scipy.stats import t


SPECTRA = ('TT', 'EE', 'BB', 'TE')


def bin_native(values):
    values = np.asarray(values, dtype=float)
    if values.shape[-1] < 1000:
        raise ValueError('Need multipoles through999')
    ell = np.arange(values.shape[-1])
    lo, hi = np.arange(25, 1000, 25), np.arange(50, 1025, 25)
    weights = 2 * ell + 1
    return np.stack([np.sum(values[..., lower:upper] * weights[lower:upper], axis=-1) /
                     np.sum(weights[lower:upper]) for lower, upper in zip(lo, hi)], axis=-1)


def gaussian_bin_variance(reference_C_ell):
    """Single-realization diagonal reference variance used for fit weights."""
    cl = np.asarray(reference_C_ell)
    ell = np.arange(cl.shape[1])
    modes = 2 * ell + 1
    factors = ell * (ell + 1.) / (2 * np.pi)
    variances = np.array([2 * cl[0]**2, 2 * cl[1]**2, 2 * cl[2]**2,
                          cl[0] * cl[1] + cl[3]**2]) / modes
    return np.array([np.sum(variances[:, lower:upper] * factors[lower:upper]**2 * modes[lower:upper]**2, axis=1) /
                     np.sum(modes[lower:upper])**2 for lower, upper in zip(range(25, 1000, 25), range(50, 1025, 25))]).T


def stable_student_statistic(residual, scatter, n):
    sem = np.asarray(scatter) / np.sqrt(n)
    # A zero scatter is not evidence for an infinitely precise stochastic pilot.
    if np.any(sem <= 0) or not np.isfinite(sem).all():
        raise ValueError('Nonpositive/nonfinite pilot scatter; investigate duplicate/degenerate realizations')
    return np.asarray(residual) / sem


def evaluate(raw_C_ell, reference_C_ell, policy):
    raw = np.asarray(raw_C_ell, dtype=float)
    reference = np.asarray(reference_C_ell, dtype=float)
    n = raw.shape[0]
    if raw.ndim != 3 or raw.shape[1:] != (4, 1025) or n != 10 or reference.shape != (4, 1025):
        raise ValueError('Expected ten4x1025 native spectra and explicit4x1025 reference')
    if not np.isfinite(raw).all() or not np.isfinite(reference).all():
        raise ValueError('Nonfinite spectrum')
    ell = np.arange(1025)
    factor = ell * (ell + 1.) / (2 * np.pi)
    raw_D_b, ref_D_b = bin_native(raw * factor), bin_native(reference * factor)
    mean, sd = raw_D_b.mean(0), raw_D_b.std(0, ddof=1)
    variance = gaussian_bin_variance(reference)
    if np.any(variance <= 0):
        raise ValueError('Invalid diagonal reference fit weights')
    fit_weights = ref_D_b / variance
    denominator = np.sum(ref_D_b * fit_weights, axis=1)
    amplitudes = np.sum(raw_D_b * fit_weights[None], axis=2) / denominator
    amplitude_mean, amplitude_sd = amplitudes.mean(0), amplitudes.std(0, ddof=1)
    amplitude_t = stable_student_statistic(amplitude_mean - 1., amplitude_sd, n)
    bin_t = stable_student_statistic(mean - ref_D_b, sd, n)
    alpha = policy['primary_familywise_alpha']
    amplitude_alpha = alpha * policy['primary_amplitude_alpha_budget_fraction']
    shape_alpha = alpha * policy['primary_shape_alpha_budget_fraction']
    amp_limit = float(t.ppf(1 - amplitude_alpha / (2 * 3), df=n - 1))
    shape_limit = float(t.ppf(1 - shape_alpha / (2 * 3 * 39), df=n - 1))
    primary_indices = [SPECTRA.index(name) for name in policy['primary_gate_spectra']]
    amp_pass = bool(np.all(np.abs(amplitude_t[primary_indices]) <= amp_limit))
    shape_pass = bool(np.all(np.abs(bin_t[primary_indices]) <= shape_limit))
    # This is stable at the TE zero-crossings and has a physical CMB cross scale.
    scale = np.sqrt(ref_D_b[0] * ref_D_b[1])
    te_residual = (mean[3] - ref_D_b[3]) / scale
    metrics = {}
    for i, name in enumerate(SPECTRA):
        entry = {'amplitude_mean': float(amplitude_mean[i]), 'amplitude_mean_minus1_pct': float(100 * (amplitude_mean[i] - 1)),
                 'amplitude_sample_SD_ddof1': float(amplitude_sd[i]), 'amplitude_mean_SEM': float(amplitude_sd[i] / np.sqrt(n)),
                 'amplitude_student_t_df9': float(amplitude_t[i]),
                 'maximum_absolute_bin_student_t_df9': float(np.max(np.abs(bin_t[i]))),
                 'bin_shape_is_formal_gate': name != 'BB'}
        if name == 'TE':
            entry['maximum_absolute_residual_over_sqrt_TT_EE'] = float(np.max(np.abs(te_residual)))
            entry['TE_ratio_used'] = False
        else:
            ratio = mean[i] / ref_D_b[i]
            entry.update({'median_bin_ratio_minus1_pct': float(100 * np.median(ratio - 1)),
                          'minimum_bin_ratio': float(ratio.min()), 'maximum_bin_ratio': float(ratio.max())})
        metrics[name] = entry
    report = {'primary_spectral_gate': 'PASS' if amp_pass and shape_pass else 'REVIEW_REQUIRED',
              'primary_amplitude_gate_pass': amp_pass, 'primary_shape_gate_pass': shape_pass,
              'primary_familywise_alpha': alpha, 'amplitude_t_limit_df9': amp_limit,
              'shape_t_limit_df9': shape_limit, 'metrics': metrics,
              'covariance_note': 'No39x39 sample covariance inversion. Empirical per-realization amplitude scatter retains bin correlation for that fitted mode. Shape t/Bonferroni uses approximate normal bandpowers.',
              'BB_note': policy['BB_note'], 'BB_is_formal_pure_CAMB_acceptance_gate': False}
    arrays = {'raw_D_b': raw_D_b, 'mean_D_b': mean, 'sample_std_D_b': sd, 'reference_D_b': ref_D_b,
              'reference_gaussian_single_realization_bin_variance': variance,
              'raw_fitted_amplitude': amplitudes, 'amplitude_student_t': amplitude_t,
              'bin_student_t': bin_t, 'TE_residual_over_sqrt_TT_EE': te_residual}
    return report, arrays


def numerical_iteration_control(main_C_ell, control_C_ell, pilot_arrays, policy):
    ell = np.arange(1025)
    factor = ell * (ell + 1.) / (2 * np.pi)
    main, control = bin_native(main_C_ell * factor), bin_native(control_C_ell * factor)
    residual = np.abs(control - main)
    scale = main.copy()
    scale[3] = np.sqrt(main[0] * main[1])
    sem = pilot_arrays['sample_std_D_b'] / np.sqrt(10)
    if np.any(scale <= 0) or np.any(sem <= 0):
        raise ValueError('Invalid numerical-control scale or pilot SEM')
    relative, sem_fraction = residual / scale, residual / sem
    max_relative, max_sem_fraction = float(relative.max()), float(sem_fraction.max())
    passed = max_relative < policy['numerical_max_fractional_auto_or_TE_correlation_scale'] and max_sem_fraction < policy['numerical_max_fraction_of_pilot_mean_SEM']
    return {'status': 'PASS' if passed else 'REVIEW_REQUIRED', 'same_native_map': True,
            'lmax': 1024, 'primary_iter': 3, 'control_iter': 1,
            'maximum_absolute_fractional_auto_or_TE_correlation_scale': max_relative,
            'maximum_fraction_of_pilot_mean_SEM': max_sem_fraction,
            'TE_denominator': 'sqrt(main_TT_D_b*main_EE_D_b), not TE_D_b',
            'per_spectrum_maximum_relative': dict(zip(SPECTRA, relative.max(1).tolist()))}
