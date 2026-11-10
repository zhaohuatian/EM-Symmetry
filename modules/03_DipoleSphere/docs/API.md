# API reference

## Shared configuration

| Field | Meaning |
|---|---|
| `lambda` | Positive vacuum wavelength [m] |
| `radius` | Sphere radius [m] |
| `epsilonSphere` | Nonzero relative epsilon, passive Im>=0 |
| `nHost` | Real positive transparent host index |
| `position` | Real Cartesian 3-vector from sphere center to source [m]; norm must exceed radius |
| `nmax` | Optional positive multipole cutoff, 1..500; refine to check convergence |

The sphere center is the origin. Translate coordinates before calling if necessary. `p` is a nonzero complex 3-vector of peak electric moment [C m], in the same global Cartesian axes. Default nmax accounts heuristically for size, source radius and surface gap. Automatic selection exceeding 500 is rejected instead of silently truncating a demanding case. Positive frequency, passive nonmagnetic homogeneous sphere and lossless exterior are required.

## `M = mlsphere.mie(C)`

Position is not required for this source-free API. If nmax is absent it is selected from sphere size. Outputs are order `l`, standard electric `a`, magnetic `b`, `Qsca`, `Qext`, `Qabs`, and nmax. Efficiencies refer to a plane wave, not to a dipole's enhancement. Its internal regular/outgoing convention is sE=-a, sM=-b, as derived in THEORY.md.

## `G = mlsphere.decay(p,C)`

Outputs `Ftotal/Frad/Fabs` are normalized to vacuum power. `FtotalHost/FradHost` are normalized to infinite homogeneous host power. Multiply host-normalized factors by nHost to obtain vacuum-normalized values. `radial/tangential` contain the two orientation-resolved host-normalized rates. `Ptotal/Prad/Pabs/P0` are in W. Gamma/Gamma0 are in s^-1 **only when p is interpreted as a transition moment**; classical peak-power conversion has the factor of four documented in the theory.

`nmax` is returned to record the calculation. A small negative Fabs in a lossless case is a numerical residual and should be checked, not interpreted as gain or hidden with a clamp.

## `F = mlsphere.farfield(p,C,theta,phi,r)`

Angles are **global** spherical angles in radians, theta in [0,pi], phi arbitrary periodic. They must have matching array shape. r defaults to 1 m and is measured from the sphere center; choose a genuine far-field distance compared with source displacement, particle size and wavelength. It evaluates the 1/r asymptote, not the near field.

`F.E` has shape `[size(theta),3]`, V/m. `F.dPdOmega` is W/sr; `F.pattern` is max-normalized for plotting. It includes direct source plus scattered field and their interference. Exact polar observation axes are handled by harmonic limits.

## `R = mlsphere.radiative(p,C,N)`

Independent full-sphere angular integration: N midpoint polar-cosine samples, 2N azimuth samples, N>=8, default 80. Returns `Prad`, `P0`, vacuum-normalized `Frad` and the sampled `Far` structure. Increase N separately from nmax; convergence in one does not establish convergence in the other.

Low-level `model`, `modes` and `angular` implement geometry/frame validation, stable source-point radial products and the m=0,+/-1 harmonics in the source-aligned frame. They are exposed for research inspection but the high-level functions are the intended entry points.
