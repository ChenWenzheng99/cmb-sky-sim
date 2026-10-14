# Independent completed-theory review

Status: PASS for artifact/convention/source/normalization checks. This is
not recovery of the original Sehgal CAMB ini or validation of map production.

Theory job225108 completed0:0 on aliws006 in51s; input-validation
job225109 completed0:0 on aliws021 in15s. Four successful artifacts are
copied to `theory/outputs`; their byte sizes and SHA256 values match the
job manifest. The failed first export225058 is retained separately and
is not treated as a successful theory result.

The actual calculation uses Python3.12.11, CAMB1.6.1 and NumPy1.26.4.
Installed model/results/correlations source hashes match all three
archived official CAMB1.6.1 sources. The executed builder hash is
`3977a4c0ddb5f92098b14d6e4ef9dc12049adc6f445e3e0e0a7627a47c72c6f2`;
the configuration file and canonical hashes retain the frozen values.

The published parameters are h=.71, Omega_b=.044, Omega_m=.264,
n_s=.96 and sigma8=.8. The actual CAMB parameter record independently
matches H0=71, omega_b=.0221804, omega_c=.110902, tau=.08,
TCMB2.7255K, YHe.24, massless N_eff3.046, pivot.05/Mpc and r0.
The latter settings remain explicit completion assumptions, not a
recovered original author ini. CAMB flatness includes radiation:
Omega_de=.7359170064922756, rather than imposing the rounded.736.
Original Smith2003 Halofit is recorded in the actual parameter string.

Two full solves with identical numerical precision yield
As=2.1500364844398247e-9 and sigma8=.8000000000000018, a relative
normalization residual2.220446049250313e-15. The independently rebuilt
As_trial*(.8/sigma8_trial)^2 exactly matches the exported amplitude.
Nonlinear lensing is recomputed at that amplitude.

All arrays are finite with the declared float64 shapes and zero-based
ell indexing. Main raw C_ell/D_ell tables are6501x4 in TT,EE,BB,TE
order and thermodynamic microkelvin squared; the DAT files roundtrip
bitwise to their NPZ arrays. Primary BB is exactly zero, autos are
nonnegative, signed TE is preserved and the T/E covariance is positive
semidefinite. Phi PP=[L(L+1)]^2 C_phi/(2pi) is dimensionless and its
conversion is exact. The theory support6500/internal6650/eta133000,
producer primary5999, phi6143 and auxiliary output1024 are distinct
and independently checked. The cutoff arrays vanish exactly above
their stated input cutoffs.

Auxiliary raw C0/D0 are exactly zero. The special unmodified library
ell0 values are preserved in NPZ and typed metadata; no positive-ell
correlation value is altered by that convention. Zero-phi D_ell is
bitwise identical to the cutoff input primary. Gauss-Legendre full-theta
sampling1.4 versus2.0 has maximum relative BB difference6.4569e-10
in the39 width25, (2ell+1)-weighted D_ell bins spanning25..999.
This measures quadrature stability for this control, not the complete
CAMB/map accuracy.

The cutoff/full-theta control versus full CAMB BB prediction differs
by median+0.0027300 percent across the same bins, with maximum absolute
fraction2.03834e-4. This comparison changes both bandlimit and numerical
lensing implementation; it does not isolate a physical cutoff cause.
No comparison with actual newly produced maps or measured mmDL PP is
asserted here. The machine-readable details and exact hashes are in
`independent_theory_review.json`.
