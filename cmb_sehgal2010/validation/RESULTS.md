# Native 2048 Sehgal-parameter CMB validation

Ten independent primary realizations (IDs 00000–00009) were remapped with
their paired stored mmDL convergence maps. Both generation and measurement
use native HEALPix NSIDE=2048. The measurement uses `lmax=1024`, three
harmonic iterations, no pixel-window correction, beam, mask or resolution
conversion. The diagnostic shows 39 bins of width 25 spanning
`25 <= ell < 1000`; each bin averages `D_ell=ell*(ell+1)*C_ell/(2*pi)` with
weights `2*ell+1`. Shading is the sample standard deviation across ten maps
(`ddof=1`), not the uncertainty on the mean.

The main reference is the self-consistent CAMB prediction with the published
Sehgal2010 parameters and the completion assumptions recorded in
`theory/sehgal2010_config.json`. A second reference applies the producer's
primary/deflection cutoffs. A third uses the same primary spectrum and the
measured mean two-point spectrum of the actual ten mmDL convergence maps in
CAMB's correlation-function lensing calculation. This conditional reference
is distinct from the self-consistent cosmological prediction. No lensing
amplitude is fitted or applied to the production maps.

| Spectrum | Fitted mean amplitude relative to self-consistent CAMB | Relative to conditional actual-kappa two-point reference |
|---|---:|---:|
| TT | +0.0471% | +0.0578% |
| EE | −0.0525% | −0.00480% |
| TE | +0.2088% | +0.0348% |
| BB | −0.6662% | −0.0395% |

These amplitudes are descriptive fits performed only in validation; each
realization's amplitude is calculated with the fixed reference weights,
then the ten amplitudes are averaged. They are not uniform averages of bin
ratios. The frozen TT/EE/TE amplitude and shape checks pass. The numerical
iteration control (same map, iterations 1 versus 3) has maximum fractional
auto-spectrum or normalized TE difference `6.3947e-8`, negligible relative
to the measured scatter.

BB is close in fractional amplitude to the self-consistent prediction, but
its remaining offset is statistically resolved by this ten-map pilot:
amplitude Student-t = −13.27 using the empirical per-map amplitude scatter.
The median BB bin difference is −0.3913%, with bin differences ranging from
−1.7195% to +1.0226%. Relative to the conditional actual-kappa reference,
the BB amplitude Student-t is −0.779 and the median bin difference is
+0.0159%. This supports a residual lens-input difference; it does not
establish equality to the pure cosmological CAMB lensing prediction or a
unique physical cause. The stored mmDL convergence is retained without
rescaling. No formal pure-CAMB BB acceptance is claimed.

The four-spectrum diagnostic is `figures/native_cmb_pilot_10.pdf` (also
PNG); its numeric data and complete reference metrics are in
`pilot_summary/native_cmb_pilot_10.npz`, `.json` and `_binned.csv`.
`validation_policy.json` records the frozen tests. Independent code/data
review is stored under `evidence/pilot_review`.

The spectral validation scope is these ten maps and the stated multipole
range. The campaign audit separately checks all 500 outputs and preservation
of the original inputs; it does not measure all 500 spectra or validate the
high-multipole end of the production maps.
