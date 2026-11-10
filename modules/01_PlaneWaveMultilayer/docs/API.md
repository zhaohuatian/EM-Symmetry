# API reference

Add the module root with `addpath`. MATLAB discovers the `+mlpw` namespace automatically.

## `R = mlpw.solve(C,method)`

| Input | Meaning |
|---|---|
| `C.lambda` | Positive vacuum wavelength, m |
| `C.epsilon` | Relative permittivity vector, left-to-right, first/last half spaces |
| `C.thickness` | Positive finite-layer thicknesses, m; length L-2; use `[]` for one interface |
| `C.kparallel` | Scalar conserved transverse wavevector, 1/m |
| `C.polarization` | `'TE'` or `'TM'` |
| `method` | Optional `'smm'` default, or `'tmm'` |

Unit left incidence, no right incidence. `R.S` is the full complex two-port matrix; `R.r=S11`, `R.t=S21`. `reflectance/transmittance/absorptance` are power fractions only for transparent propagating ports. `R.q`, `R.z`, `R.Aleft`, `R.Bright` support field reconstruction. `R.S_smm` always records the stable cascade; `R.transfer` exists for TMM. Internal fields use stable S coefficients even in TMM mode, explicitly separating a TMM global comparison from robust field evaluation.

## `F = mlpw.field(R,z)`

`z` is a vector of global positions in m, with first interface at z=0. The final interface is at sum(C.thickness). `F.E` and `F.H` are N-by-3 Cartesian arrays. For an incident amplitude of 1 V/m they are in V/m and A/m; `F.Sz` is W/m². Scale the returned field linearly for another incident electric amplitude. An interface sample belongs to the layer on its right; use one-sided points to check continuity. Tangential E/H are continuous, while normal E need not be.

`mlpw.interface` and `mlpw.cascade` are low-level functions documented in the derivation; users usually need only `solve` and `field`. Singular exact light-line points are rejected; sample their limits rather than inserting a zero denominator.
