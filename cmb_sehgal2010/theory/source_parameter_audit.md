# Primary-source parameter and theory audit

This audit supports the independent `sehgal2010_config.json` calculation.
It preserves archived maps and does not establish recovery of the original
author CAMB parameter file.

| Item | Value or conclusion | Primary evidence |
|---|---|---|
| Published cosmology | Omega_b=.044, Omega_m=.264, rounded Omega_Lambda=.736, h=.71, n_s=.96, sigma8=.80 | [Sehgal et al.2010, Section1, printedpage2](https://arxiv.org/abs/0908.0540). Snapshot `source_evidence/Sehgal2010_simulations.pdf`. |
| Physical densities used | omega_b=.0221804, omega_c=.110902, H0=71 | Algebra from the published tuple with the explicitly chosen massless-neutrino completion. |
| Approximate gas optical depth | tau_T approximately.08 for the specified z=10 ionization scenario | [Same paper, Section2.4.2](https://arxiv.org/pdf/0908.0540). This motivates a completion; it does not identify the primary-CMB CAMB tau. |
| Author primary spectrum | Headerless `bode_almost_wmap5_lmax_1e4_scalCls.dat`; no recoverable full CAMB ini in the archived loader | [Pinned author resource](https://github.com/dwhan89/cosmikyu/blob/927c56f76a9c9f7c517a86a52586b8958c308328/resources/bode_almost_wmap5_lmax_1e4_scalCls.dat), [pinned loader](https://github.com/dwhan89/cosmikyu/blob/927c56f76a9c9f7c517a86a52586b8958c308328/cosmikyu/sehgal.py). |
| Author TCMB constant | 2.726K is used in foreground SED/unit conversion | Same pinned `sehgal.py`; it is not evidence for the CAMB background TCMB value. |
| Unrecovered primary settings | Exact tau, neutrinos, N_eff, YHe, scalar pivot and original numerical ini | No such settings in the retained headerless table/loader. No exact-reproduction claim is made. |
| Explicit completion | tau.08, mnu0, massive-neutrino count0, N_eff3.046, YHe.24, TCMB2.7255K, pivot.05/Mpc, r0, no running, flat LambdaCDM | Calculation choices recorded individually in the JSON; these are not attributed to the authors. |
| Nonlinear phi reference | Original Smith2003 Halofit in CAMB | [Smith et al.2003](https://arxiv.org/abs/astro-ph/0207664), [CAMB1.6.1 nonlinear model](https://github.com/cmbant/CAMB/blob/1.6.1/camb/nonlinear.py). This is an explicit continuum approximation, not identity with a finite-box/GAN map. |
| Scalar normalization | Derive As from linear sigma8(z0)=.8 using two complete solves at identical final numerical precision | [CAMB1.6.1 result API](https://github.com/cmbant/CAMB/blob/1.6.1/camb/results.py). The final solve recalculates nonlinear lensing; no stale nonlinear transfer rescaling. |
| Raw CMB arrays | Thermodynamic microkelvin squared; zero-based ell; columns TT, EE, BB, TE | Same CAMB1.6.1 result getter source. Scalar-only unlensed BB is zero. |
| Phi normalization | PP=[L(L+1)]^2 C_L(phi,phi)/(2pi), dimensionless | [CAMB1.6.1 result source](https://github.com/cmbant/CAMB/blob/1.6.1/camb/results.py), [correlation function source](https://github.com/cmbant/CAMB/blob/1.6.1/camb/correlations.py). No TCMB-squared factor enters phi. |
| Auxiliary continuum relensing | Full theta Gauss-Legendre integral; primary5999 and PP6143; output1024 | Same version-pinned correlation function source; uses no CAMB2.0-only keywords. |

The original paper's rounded Omega_Lambda is not imposed in addition to
flatness: CAMB includes radiation in its background sum. Actual z=0
densities are exported rather than asserting the rounded tuple as an
exact background equality.

The archived prior comparison reached sigma8=.8000808913306604 with
As=2.1504713402584196e-9. That result is a separate historical control;
this builder recomputes As and enforces its declared1e-6 relative
normalization tolerance at unchanged final numerical settings.

There are three distinct references: the full self-consistent completed
Sehgal cosmology, its producer-bandlimit CAMB-phi control, and the same
bandlimited primary evaluated with the measured mmDL mean PP. The last
conditions on the measured input power; it does not imply a recovered
original cosmology, independent theory validation, or a unique explanation
of any BB residual. The PP-only Gaussian-deflection reference omits
connected higher-order statistics of non-Gaussian convergence.

Source filenames, URLs and SHA256 values are recorded in
`source_evidence/source_manifest.json`. No large CAMB run or map
production is asserted by this source audit.
