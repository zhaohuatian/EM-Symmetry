# API and numerical controls

Dependency: module 01_PlaneWaveMultilayer must be on the MATLAB path. From the collection root run `setup`. The rt function delegates the planar stack calculation to mlpw.solve; it only maps source geometry and reference phases.

## Shared inputs

`p` is a nonzero complex 3-vector of peak electric dipole moment in C m. `D` uses the original author's geometry fields:

| Field | Definition |
|---|---|
| `WL0` | Vacuum wavelength [m] |
| `Eplist` | Relative epsilon vector, **bottom to top** |
| `dl` | Strictly increasing global interface coordinates [m], length L-1 |
| `nDn` | Source layer minus one; integer 0..L-1 |
| `POSD` | Source z coordinate [m], strictly inside its layer |

The source host must have real positive epsilon. Passive absorbing other layers are allowed. For a homogeneous half-space with no interfaces use a one-element Eplist, empty dl and nDn=0. The present release uses its own stable RT core and has no dependency on the collection's plane-wave folder.

## `F = mldipole.farfield(p,D,theta,phi,side,r)`

Theta and phi have equal array shape, in radians. Theta is 0<=theta<pi/2 measured from the **outward axis of the chosen side**. Phi uses the global x/y axes. Side defaults to +1 (upper); -1 is lower. r defaults to 1 m and must be an appropriate far-field distance for a physical interpretation.

`F.E` has size `[size(theta),3]` and contains Cartesian complex E in V/m. `F.dPdOmega` is W/sr, independent of chosen r. `F.pattern` is max-normalized for plotting only. Exterior observation epsilon must be real positive. Exact host/interior critical samples should be approached by limits rather than forcing zero denominators.

## `R = mldipole.radiative(p,D,N)`

Two hemispheres, N midpoint samples in cos(theta) and 2N periodic azimuth samples. N defaults to 160, integer >=8. Outputs: `PUp`, `PDn`, `Prad`, `P0` in W, `Frad=Prad/P0`, and `Up/Dn` far-field grids. This is an escaped-power definition. Increase N and inspect critical-angle convergence; the default is not a universal error bound.

## `G = mldipole.decay(p,D,opt)`

Optional `opt.eta` default .02, `opt.RelTol` default 1e-7, `opt.AbsTol` default 1e-9 and `opt.umax` (dimensionless transverse-wavevector cutoff, automatically chosen from source-interface distance if absent). It integrates the scattered self-Green tensor and adds the analytic uniform-host imaginary Green term.

Outputs: `Ftotal` vacuum-normalized; `Ptotal/P0` in W; `Gsc` in 1/m; `ImG` is the regularized imaginary tensor; `Esc` is the scattered feedback E in V/m. `Gamma/Gamma0` are in s^-1 only if the numerical input p is interpreted as a **quantum transition** dipole. Otherwise use the classical power outputs.

## `pat = mldipole.pattern_grid(p,D,objective,side)`

Compatibility entry point for original `kx_grid/ky_grid` arrays. `objective=0` yields W/sr, `objective=1` yields dP/d(ux)d(uy) with ux=kx/kout. Samples outside the exterior light cone are NaN. No optical objective transmission or camera response is modeled.

Low-level `mldipole.model` and `mldipole.rt` implement validation and source-plane stack responses. End users generally call the three high-level field/power functions.
# Diagonal-anisotropic far-field extension

In addition to the isotropic `D` interface documented below, a separate `C` interface uses `lambda` in meters, `epsilonDiagonal` as an L-by-3 bottom-to-top array, `interfaces` as L-1 increasing global z positions, and `sourceZ` in meters. Both exteriors must be isotropic; the observation exterior must be transparent. Call `F=mldipole.reciprocal_farfield(p,C,theta,phi,side,r)` with side=+1/-1 and theta from that side's outward normal. It returns E [V/m], dPdOmega [W/sr], pattern and residual diagnostics. This is far field only.

`S=mlpw.diagonal_stack(C,kx,ky,Ain,Bin)` solves auxiliary plane-wave fields; kx and ky are dimensional wavevectors [1/m], and incoming electric coefficients Ain/Bin are 2-by-K arrays at bottom/top interface references. Each column is a separate illumination, with rows TE/TM. `S.coefficients` is 4-by-L-by-K: two upward modes followed by two downward modes. `S.q` is normalized to k0. `S.lowerReference` and `S.upperReference` store their distinct origins. `mlpw.modal_field(S,z)` reconstructs auxiliary E and H at requested global z samples. No dipole-source jump is included in this auxiliary solve.

See the [walkthrough](FIGURE2_RECONSTRUCTION.md) for eigenmode definitions, light-line treatment, material provenance and full example commands.

