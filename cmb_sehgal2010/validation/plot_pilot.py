#!/usr/bin/env python3
"""Four-spectrum native2048 ten-map diagnostic; sample-SD shading only."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from validate_inputs import check_record, sha256


def plot(summary, out):
    summary, out = Path(summary), Path(out)
    metadata = json.loads(summary.with_suffix('.json').read_text())
    check_record(summary, metadata['artifact'])
    if metadata['policy']['nside'] != 2048 or metadata['policy']['ensemble_SD_ddof'] != 1:
        raise ValueError('Expected native2048 ten-map sample-SD summary')
    with np.load(summary, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    lo = np.arange(25, 1000, 25)
    ell = np.arange(1025)
    centers = np.array([np.average(ell[l:l + 25], weights=2 * ell[l:l + 25] + 1) for l in lo])
    main = 'self_consistent_CAMB'
    mean, sd = arrays[main + '__mean_D_b'], arrays[main + '__sample_std_D_b']
    ref = arrays[main + '__reference_D_b']
    if mean.shape != (4, 39) or sd.shape != (4, 39) or not np.isfinite(mean).all() or not np.isfinite(sd).all():
        raise ValueError('Invalid pilot bands')
    references = [(main, 'CAMB: Sehgal parameter completion', 'black', '--'),
                  ('producer_bandlimit_CAMB', 'CAMB with producer cutoffs', '0.5', ':')]
    matched = 'conditional_actual_mean_kappa_2point'
    if matched + '__reference_D_b' in arrays:
        references.append((matched, 'Conditional measured-κ two-point reference', '#c77718', '-.'))
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 6.4), sharex='col',
                             gridspec_kw={'height_ratios': [3, 1.5], 'hspace': .10, 'wspace': .34})
    spectrum_names = ['TT', 'EE', 'BB', 'TE']
    bounds = {}
    for index, name in enumerate(spectrum_names):
        top, bottom = axes[:, index]
        top.set_title(name)
        top.plot(centers, mean[index], color='#2465a6', lw=1.7, label='Native2048 mean of10')
        low, high = mean[index] - sd[index], mean[index] + sd[index]
        if index < 3:
            if np.any(high <= 0):
                raise ValueError('Nonpositive auto-spectrum upper SD bound')
            top.set_yscale('log')
            display_low = np.maximum(low, np.finfo(float).tiny)
        else:
            display_low = low
            top.axhline(0, color='0.8', lw=.8)
        top.fill_between(centers, display_low, high, color='#2465a6', alpha=.20,
                         label='± sample SD (ddof=1)')
        for key, label, color, style in references:
            top.plot(centers, arrays[key + '__reference_D_b'][index], color=color,
                     ls=style, lw=1.3, label=label)
        if index == 3:
            scale = np.sqrt(ref[0] * ref[1])
            residual, scatter = (mean[3] - ref[3]) / scale, sd[3] / scale
            bottom.set_ylabel(r'$\Delta D_b^{TE}/\sqrt{D_b^{TT}D_b^{EE}}$', fontsize=9)
        else:
            residual, scatter = 100 * (mean[index] / ref[index] - 1), 100 * sd[index] / ref[index]
            bottom.set_ylabel('Mean / CAMB − 1 [%]', fontsize=9)
        bottom.axhline(0, color='black', ls='--', lw=.8)
        bottom.plot(centers, residual, color='#2465a6', lw=1.5)
        bottom.fill_between(centers, residual - scatter, residual + scatter, color='#2465a6', alpha=.20)
        for key, _, color, style in references[1:]:
            alternative = arrays[key + '__reference_D_b'][index]
            delta = (alternative - ref[index]) / scale if index == 3 else 100 * (alternative / ref[index] - 1)
            bottom.plot(centers, delta, color=color, ls=style, lw=1.1)
        top.set_ylabel(r'$D_b\;[\mu\mathrm{K}_{\rm CMB}^2]$')
        bottom.set_xlabel(r'$\ell$')
        top.set_xlim(0, 1000)
        for axis in (top, bottom):
            axis.grid(alpha=.18)
            axis.tick_params(labelsize=9)
        bounds[name] = {'sample_SD_lower_minimum_uK2': float(low.min()),
                        'negative_lower_log_display_bins': int(np.count_nonzero(low <= 0)) if index < 3 else 0}
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=9, bbox_to_anchor=(.5, .012))
    fig.subplots_adjust(bottom=.18, top=.92, left=.06, right=.98)
    fig.suptitle('Sehgal2010-completion CMB pilot: native NSIDE2048, IDs0–9', fontsize=13)
    out.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ['pdf', 'png']:
        path = out / ('native_cmb_pilot_10.' + extension)
        fig.savefig(path, dpi=180)
        outputs[extension] = {'path': str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': sha256(path)}
    plt.close(fig)
    result = {'summary': {'path': str(summary.resolve()), 'sha256': sha256(summary)},
              'summary_metadata_sha256': sha256(summary.with_suffix('.json')), 'outputs': outputs,
              'script_sha256': sha256(__file__), 'source_nside': 2048, 'IDs': list(range(10)),
              'spectrum_order': spectrum_names, 'ell_display_interval': '[25,1000)',
              'bin_quantity': '2ell+1 weighted D_ell', 'shade': 'sample SD ddof1; not SEM',
              'main_reference': main, 'references': [x[0] for x in references],
              'TE_lower_panel': 'signed difference divided by sqrt(main-reference TT*EE), no TE ratio',
              'no_fitted_lensing_amplitude': True, 'shade_bounds': bounds}
    (out / 'native_cmb_pilot_10_selection.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(plot(args.summary, args.out), indent=2, allow_nan=False))
