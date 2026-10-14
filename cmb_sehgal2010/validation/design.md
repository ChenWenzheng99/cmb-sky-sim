# Independent native-grid pilot validation

This directory is separate from the generator and theory builder. The pilot
contains IDs 0–9 on the original NSIDE2048 grid. No NSIDE512 result or single-ID
curve substitutes for that ensemble. The policy is fixed before its maps and
spectra are produced; no input spectrum or lensing amplitude is fitted to make
the maps agree.

## Source and transformation audit

The original producer snapshot is
`production/original_snapshot/lensed_CMB_generation.py`, SHA256
`ab7c9bb4dfc0f1bcdfbbb04cf7658a8feb73837193618a8013067111de141d21`.
Its primary bandlimit is 5999 and its native κ transform and gradient
deflection bandlimit are 6143. The original `synalm` call uses healpy's legacy
TT,TE,EE,BB order. The release makes `new=True` and TT,EE,BB,TE explicit,
preserving the covariance while adding a reproducible, recorded task seed.
The old unrecorded seeds are not recovered.

For harmonic index L≥1 the producer uses

\[
\phi_{LM}=2\kappa_{LM}/[L(L+1)],\qquad
d_{LM}=\sqrt{L(L+1)}\phi_{LM}.
\]

It sets the φ and deflection monopoles to zero and retains the dipole.
Consequently raw auto/cross spectra must satisfy
\(C_L^{\phi\phi}=4C_L^{\kappa\kappa}/[L(L+1)]^2\),
\(C_L^{\kappa\phi}=2C_L^{\kappa\kappa}/[L(L+1)]\), and their
gradient-deflection counterparts with the same positive transfer. The
dimensionless CAMB correlation-function input is
\(PP_L=[L(L+1)]^2C_L^{\phi\phi}/(2\pi)=2C_L^{\kappa\kappa}/\pi\)
for L≥1. It is not a temperature spectrum.

`check_producer_input_spectra` verifies these relations independently from the
saved producer alms' small spectra. This adds no κ spherical transform. The
saved `deflection_transfer` means κ→d, whereas the producer's named gradient
step is φ→d. The positive auto/cross checks detect a sign reversal as well as
an amplitude or unit error.

The remapping passes three T/E/B primary alms to lenspyx, applies its spin-2
polarization treatment, and returns T/Q/U. The producer uses
`epsilon=1e-6`, the full native HEALPix geometry, no beam, no mask, no pixel
window correction, and no map downgrade. The source κ pixel directions are
retained. A missing source COORDSYS is recorded as the owner's declared G
frame; the upstream ICRS provenance conflict remains explicit. Full-sky
isotropic spectra cannot physically certify that declaration.

Official interface references: [healpy synalm ordering](https://healpy.readthedocs.io/en/stable/generated/healpy.sphtfunc.synalm.html),
[healpy map2alm](https://healpy.readthedocs.io/en/1.16.2/generated/healpy.sphtfunc.map2alm.html),
[lenspyx remapping source](https://lenspyx.readthedocs.io/en/stable/_modules/lenspyx/lensing.html),
and [CAMB correlation-function conventions](https://camb.readthedocs.io/en/latest/correlations.html).
The numerical run is pinned to CAMB1.6.1, lenspyx2.0.0 and ducc0.38;
version-specific calls are checked against the pinned source, not assumed
from the current documentation.

## Theory inputs and references

Raw primary/lensed Cℓ columns are explicitly TT,EE,BB,TE in thermodynamic
µK_CMB². Dℓ is obtained by multiplying ℓ(ℓ+1)/(2π) exactly once. The complete
6500 table and the producer primary cut at5999 are separate arrays. Scalar
unlensed BB is exactly zero and r=0. Published Ωb=.044, Ωm=.264, h=.71,
ns=.96 and σ8(z=0)=.8 follow [Sehgal et al.](https://arxiv.org/abs/0908.0540).
The missing original CAMB τ, neutrino, helium, TCMB and recombination settings
remain declared completion assumptions; the default mmDL SED temperature
2.726K does not recover the original CAMB TCMB input.

The main reference is self-consistent CAMB lensed CMB for the declared
completion. Its σ8 normalization uses two full solves with the same numerical
settings and recomputes nonlinear lensing at the calibrated As. A separate
continuum control cuts the same primary at5999 and CAMB PP at6143, using the
full angular range. Sampling factors1.4 and2 and zero-φ identity are recorded.
These controls isolate integration and retained-bandlimit effects.

After the ten maps exist, `summarize_kappa_inputs.py` extracts their actual
mean PP from the producer's existing alms. The separate theory helper can
relens the same primary with this PP. This is a conditional two-point
Gaussian-deflection reference, not an independent cosmological prediction
or the mean of ten individually relensed predictions. It does not include
all connected higher-order statistics of non-Gaussian mmDL κ.

No narrow κ multipole power ratio alone determines a BB discrepancy: lensing
BB couples many E and lensing modes. A control changing only L<25 PP can
isolate the effect of the low-L deficit if needed. An earlier three-κ/FFP
primary calculation gives a prior matched-reference median BB difference
−0.55233% relative to Sehgal CAMB PP, and maximum0.006085% change when PP
above4096 is removed. Those values, saved in `prior_BB_expectation.json`, are
conditional historical estimates, not a forecast certified for the new ten
maps and their new primary.

## Measurement and predeclared statistics

Each native T/Q/U map is measured with joint spin transforms at lmax1024,
iter3, uniform quadrature (`use_weights=False`, `use_pixel_weights=False`).
TT,EE,BB,TE are extracted from the T/E/B alms. Display bins are [25,50),
[50,75), …, [975,1000), weighted by2ℓ+1 after Dℓ conversion. Their x positions
are the same weighted mean ℓ. No beam, pixel window, monopole/dipole removal,
mask, or coordinate change is applied during measurement. Full source
hashes, FITS contracts and all150,994,944 finite pixel values are checked.

Plots show the mean of ten maps and sample SD with ddof=1. SEM=SD/√10 is used
only in explicitly named mean-consistency tests. Ten realizations cannot
provide an invertible unrestricted39×39 sample covariance.

For TT,EE,TE the primary familywise α=.01 is split equally between three
per-realization fitted amplitude tests and117 bin shape tests. Amplitude fit
weights are precomputed from diagonal Gaussian reference variance, while
the fitted amplitudes' empirical scatter preserves bin correlations for
that mode. Student-t df9 thresholds and Bonferroni are approximate bandpower
diagnostics, not exact non-Gaussian likelihood claims. Signed TE residuals
are divided by √(reference TT×reference EE); no TE zero-crossing ratio is used.

BB has an amplitude, median/min/max ratio and empirical significance report
against every named reference. It is not forced through the primary gate.
Good overall agreement requires explicit inspection of its actual residuals
and the matched-κ diagnostic. This separation matters because lensing BB has
correlated, non-Gaussian sample variance; see
[Smith, Hu and Kaplinghat](https://arxiv.org/abs/astro-ph/0607315).

ID0 is also transformed at the same1024 cutoff with iter1. Predeclared
iteration thresholds are maximum absolute auto fractional change (or TE
correlation-scale change)<1e-4 and maximum change<0.01 of the pilot-mean SEM.
The previous native-grid BB iter1/3 control differed by only5.8e-8 in
fractional binned power. The anticipated Gaussian mean-SEM fraction near
the upper bin is about.002, so the SEM threshold corresponds to roughly2e-5
fractional power, over300 times that old numerical difference. These are
numerical controls, not tolerances to hide a lens-model discrepancy.

## Commands

Run from the release root. The root agent schedules native generation and
measurement; this validation directory itself submits no jobs.

```bash
python -m unittest discover -s validation -p test_validation.py -v
python validation/validate_inputs.py --root . --report validation/results/theory_input_audit.json
python validation/measure_native_cmb.py --id 0 --map /path/00000/lensedcmb_new_00000.fits --out /path/pilot_spectra
python validation/summarize_kappa_inputs.py --producer-dir /path/mmDL_maps_2048 --theory theory/outputs/sehgal2010_theory.npz --out validation/results/kappa_inputs
python theory/relens_primary_with_pp.py --theory theory/outputs/sehgal2010_theory.npz --pp-npz validation/results/kappa_inputs/producer_kappa_ensemble_10.npz --out validation/results/matched_reference
python validation/summarize_pilot.py --root . --spectra-dir /path/pilot_spectra --matched-reference validation/results/matched_reference/matched_mmdl_pp_reference.npz --out validation/results/pilot
python validation/plot_pilot.py --summary validation/results/pilot/native_cmb_pilot_10.npz --out validation/results/pilot/figures
```

The first nine small analytic-array tests validate conventions and rejection
paths only. They neither generate a lower-resolution sky nor certify the
native remapping. Native numerical status and scientific BB interpretation
are reported separately when the actual pilot completes.
