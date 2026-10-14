# Sehgal2010 primary-CMB theory

This independent directory defines the CAMB input for a replacement primary
CMB component. Old maps and the prior data-paper release are preserved.

## Published parameters and missing settings

[Sehgal et al. (2010), Section 1](https://arxiv.org/abs/0908.0540) gives
(Omega_b, Omega_m, Omega_Lambda, h, n_s, sigma_8) =
(0.044, 0.264, 0.736, **0.71**, 0.96, 0.80). Thus H0=71, omega_b=.0221804
and omega_c=.110902 for massless neutrinos. The paper's gas model assumes
hydrogen reionization at z=10 and approximately tau=.08; that does not
recover the primary-CMB CAMB optical-depth setting.

The pinned author resource
[bode_almost_wmap5_lmax_1e4_scalCls.dat](https://github.com/dwhan89/cosmikyu/blob/927c56f76a9c9f7c517a86a52586b8958c308328/resources/bode_almost_wmap5_lmax_1e4_scalCls.dat)
has no parameter header. Its associated
[author code](https://github.com/dwhan89/cosmikyu/blob/927c56f76a9c9f7c517a86a52586b8958c308328/cosmikyu/sehgal.py)
loads those spectra but supplies no complete original CAMB ini. The
code's 2.726 K SED temperature does not establish the CAMB background
temperature. Exact original tau, neutrino settings, YHe and scalar pivot
remain unrecovered.

The explicit completion in `sehgal2010_config.json` is tau=.08, massless
neutrinos, N_eff=3.046, YHe=.24, TCMB=2.7255 K, pivot=.05/Mpc, r=0,
zero running and flat LambdaCDM. These values define the generated
calculation; they are not described as a recovered original parameter
file. With r=0 the scalar primary BB is exactly zero. Lensing generates
BB from the remapped polarization.

CAMB includes radiation when enforcing flatness, so the calculated
Omega_Lambda is slightly below the paper's rounded .736; the actual
background densities are exported. Original Smith2003 nonlinear
lensing is specified as an explicit continuum-reference choice.

## Builder and array contract

Validate the configuration without running CAMB:

```bash
python build_sehgal2010_theory.py --print-config
```

The high-L calculation must be run in its coordinated compute allocation,
with the configured CAMB version:

```bash
python build_sehgal2010_theory.py --config sehgal2010_config.json --out outputs
```

The builder runs a normalization trial and the final calculation with
identical numerical settings. It derives
As_final=As_trial*(sigma8_target/sigma8_trial)^2 and recomputes nonlinear
lensing completely at the final amplitude. It does not fit the archived
BB or reuse nonlinear transfers after changing As. The achieved linear
sigma8 is checked against .8 with relative tolerance 1e-6.

Output `outputs/sehgal2010_theory.npz` has this explicit zero-based contract:

| Key | Shape | Meaning |
|---|---|---|
| `ell` | (6501,) | int64 ell=0,...,6500. |
| `unlensed_raw_cl` | (6501,4) | C_ell in thermodynamic microkelvin squared, columns TT, EE, BB, TE. |
| `lensed_raw_cl` | (6501,4) | Self-consistent CAMB scalar+lensing C_ell, same units and order. |
| `unlensed_D_ell`, `lensed_D_ell` | (6501,4) | Corresponding ell(ell+1)C_ell/(2pi). |
| `producer_unlensed_raw_cl` | (6501,4) | Primary raw C_ell through ell5999; exact zero above5999. |
| `phi_raw_cl` | (6501,) | Dimensionless raw C_L(phi,phi). |
| `phi_PP` | (6501,) | Dimensionless [L(L+1)]^2 C_L(phi,phi)/(2pi). |
| `producer_cutoff_ell` | (1025,) | int64 ell=0,...,1024 for auxiliary cutoff controls. |
| `producer_cutoff_lensed_raw_cl`, `producer_cutoff_lensed_D_ell` | (1025,4) | Continuum reference using primary through5999 and CAMB phi PP through6143. |
| `producer_cutoff_phi_PP` | (6501,) | CAMB phi PP with exact zero above6143. |
| `producer_cutoff_sampling2_raw_cl`, `producer_cutoff_sampling2_D_ell` | (1025,4) | Same cutoff control at quadrature sampling factor2.0 instead of1.4. |
| `producer_cutoff_zero_phi_D_ell` | (1025,4) | Zero-lensing identity control; exact equality to the input primary D_ell is checked. |
| `raw_correlation_zero_mode` | (3,4) | Unmodified special ell0 return from the cutoff, sampling2.0 and zero-phi correlation calls; not physical D_ell. |

These conventions follow the
[CAMB 1.6.1 result getters](https://github.com/cmbant/CAMB/blob/1.6.1/camb/results.py).
The standalone raw-Cl DAT files are explicitly labeled; they should not be
passed to a reader that assumes CAMB D_ell tables. `theory_manifest.json`
records the full configuration, calibrated As, achieved sigma8, array
checks, exact actual CAMB parameters, versions and file/source hashes.

The producer preserves primary lmax5999 and kappa/deflection lmax6143 at
NSIDE2048. Table lmax6500 provides theory support above both cutoffs;
it does not raise the map producer's primary bandlimit. CAMB's separate
`lens_margin=150` makes its internal max_l6650. Those meanings follow
the [pinned set_for_lmax implementation](https://github.com/cmbant/CAMB/blob/1.6.1/camb/model.py).
Numerical boosts2, lens-potential accuracy2, matter kmax5/Mpc and
max_eta_k=20*6650=133000 are recorded choices. These are not a standalone
convergence certification.

The new map generation must draw `unlensed_raw_cl` with an explicit
primary cutoff5999 and healpy `new=True` order TT, EE, BB, TE, then use
the selected existing mmDL kappa map. CAMB `phi_raw_cl` is only the
self-consistent continuum reference, not a replacement map-lensing
input. Consequently the new primary cosmology changes TT/TE/EE as well
as lensing BB; it is distinct from the preserved FFP10-primary maps.

The auxiliary producer-cutoff arrays use
[CAMB 1.6.1 correlations.lensed_cls](https://github.com/cmbant/CAMB/blob/1.6.1/camb/correlations.py)
with its full angular range (`theta_max=None`), Gauss-Legendre integration,
and explicitly truncated primary and CAMB phi PP. They preserve the two
producer bandlimits while leaving the main full CAMB `lensed_raw_cl`
unchanged. Sampling factors1.4 and2.0 and the zero-PP exact identity are
stored as separate numerical controls. The quadrature comparison is not
a certification of other CAMB accuracy settings or the map generator.

The pinned correlation implementation sets its internal `lfacs[0]=1`.
Its returned ell0 has a special normalization and cannot be converted
with ell(ell+1)/(2pi). The wrapper preserves that unmodified return in
the NPZ and typed manifest metadata, defines auxiliary raw C0 and D0
as zero, and leaves every positive-ell correlation result unchanged.
The main CAMB getter arrays are untouched. The same convention applies
to the measured-PP helper. A cheap check requires no CAMB calculation:

```bash
python -m unittest discover -s . -p test_theory_conventions.py
```

After the producer has supplied its measured ensemble mean phi PP, a
separate conditional reference can be built without another kappa map
transform:

```bash
python relens_primary_with_pp.py \
  --theory outputs/sehgal2010_theory.npz \
  --pp-npz /path/to/producer_kappa_ensemble.npz \
  --pp-key mean_phi_PP --ell-key ell --out matched_reference
```

The explicitly selected PP must be a float64 zero-based vector through
L6143 with normalization [L(L+1)]^2 C_L(phi,phi)/(2pi). The helper verifies
the adjacent theory manifest and input hashes, uses the same primary
cutoff and full-theta integration, and writes raw C_ell and D_ell through
1024 with separate provenance. This is a Gaussian-deflection continuum
reference conditional on that measured mean PP; it is not an independent
cosmological prediction or a fit. Higher connected statistics of the
non-Gaussian mmDL convergence are not included by the PP-only reference.

## Current status

Code/configuration and primary-source evidence are prepared. No high-L
theory artifacts are labeled complete before a successful coordinated
run and export checks. The failed first export's diagnostic log is retained
in `source_evidence/theory_225058_failed.log`; it produced no successful
theory artifact. The previous control's As=2.1504713403e-9 and achieved
sigma8=.8000808913 are not copied as this calculation's normalization.
