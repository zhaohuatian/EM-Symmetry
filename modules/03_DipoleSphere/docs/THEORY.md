# Derivation: an exterior dipole and a dielectric/metal sphere

## 1. Physical model and what is being solved

A sphere of radius a and relative permittivity epsilon_s is centered at the origin. Its host has real positive index nh, epsilon_h=nh², and permeability mu0. A maintained electric dipole p is located at r0 with |r0|>a. Both dielectric and passive metallic epsilon_s are allowed under exp(-i omega t), so Im(epsilon_s)>=0.

Use k=nh omega/c, x=ka, x0=k|r0|, y=k0 sqrt(epsilon_s) a and m=sqrt(epsilon_s)/nh. The boundary problem is exactly separable in vector spherical waves. It has electric (TM-type) and magnetic (TE-type) multipole scattering coefficients rather than a single plane-wave r/t. A scalar “sphere reflection coefficient” cannot describe all angles and polarizations.

Rotate the local z axis to r0/|r0|. The code uses a real orthogonal frame W: local vectors transform to global by W. Then p_local=W^T p, and only m_azimuth=0,+1,-1 source components are needed. This m_azimuth is distinct from the refractive-index ratio m.

## 2. Radial functions and vector spherical harmonics

Let j_l and h_l=h_l^(1) be spherical Bessel and outgoing Hankel functions. Define psi_l(z)=z j_l(z), xi_l(z)=z h_l(z). Primes below mean differentiation of the **whole argument**, not differentiation of the order or physical radius.

Choose scalar harmonics normalized by integral |Y_lm|² dOmega=1 with the Condon-Shortley phase. Set

$$\Psi_{lm}=\frac{\nabla_\Omega Y_{lm}}{\sqrt{l(l+1)}},\qquad
\Phi_{lm}=\hat r\times\Psi_{lm}.$$

For radial function z_l(kr), define

$$M_{lm}=z_l(kr)\Phi_{lm},$$
$$N_{lm}=\frac1k\nabla\times M_{lm}
=-\left[\frac{\sqrt{l(l+1)}z_l(kr)}{kr}Y_{lm}\hat r+
\frac{[kr z_l(kr)]'}{kr}\Psi_{lm}\right].$$

M has tangential electric TE character; N has electric TM character and a radial component. The angular bases are orthonormal on the sphere. On a complex radial argument the source-side Green expansion conjugates the **angular** harmonic only, not the outgoing radial function. Taking an ordinary absolute square of h in a feedback Green tensor would give the wrong total-decay expression.

## 3. Boundary matching and Mie coefficients

For every (l,m), write the incident regular field, outgoing scattered field and internal regular field. For an M-type excitation this is `M(j)+sM M(h)` outside and cM M(j_internal) inside. For N-type use `N(j)+sE N(h)` and cE N(j_internal). Enforce tangential E and H at r=a independently in each channel; no angular modes mix for a homogeneous isotropic sphere.

Introduce logarithmic derivatives

$$J_x=\psi_l'(x)/\psi_l(x),\quad H_x=\xi_l'(x)/\xi_l(x),\quad J_y=\psi_l'(y)/\psi_l(y).$$

For M waves, tangential E is z_l and tangential H is proportional to n times `[z z_l]'/z`. Their two boundary equations give

$$s_M=\frac{j_l(x)}{h_l(x)}K_M,\qquad
K_M=\frac{n_sJ_y-n_hJ_x}{n_hH_x-n_sJ_y}.$$

For N waves the roles of the radial derivative and function interchange, yielding

$$s_E=\frac{j_l(x)}{h_l(x)}K_E,\qquad
K_E=\frac{n_sJ_x-n_hJ_y}{n_hJ_y-n_sH_x}.$$

These outgoing-to-regular ratios use **sE=-a_l, sM=-b_l**, where a_l and b_l are the usual electric/magnetic Mie coefficients. Equivalently,

$$a_l=\frac{m\psi_l(y)\psi_l'(x)-\psi_l(x)\psi_l'(y)}
{m\psi_l(y)\xi_l'(x)-\xi_l(x)\psi_l'(y)},$$
$$b_l=\frac{\psi_l(y)\psi_l'(x)-m\psi_l(x)\psi_l'(y)}
{\psi_l(y)\xi_l'(x)-m\xi_l(x)\psi_l'(y)}.$$

The extra minus sign arises because the chosen spherical regular/outgoing wave coefficients are defined directly by the tangential boundary equations. It is retained consistently in all source and far-field formulas; it is not corrected by taking absolute values.

`mlsphere.mie` returns the standard a_l,b_l and plane-wave efficiencies

$$Q_{sca}=\frac2{x^2}\sum_l(2l+1)(|a_l|^2+|b_l|^2),\qquad
Q_{ext}=\frac2{x^2}\sum_l(2l+1)\operatorname{Re}(a_l+b_l).$$

For passive lossless spheres Qext=Qsca. These are plane-wave efficiencies, distinct from the power of a nearby dipole source.

## 4. Scattered Green tensor and the dipole source

The electric field of a dipole is E=mu0 omega² G p. The outgoing homogeneous Green expansion has the normalization ik times products of normalized vector spherical waves. Scattering converts a regular incident spherical wave to sM/sE times its outgoing counterpart, so for both r and r0 outside the sphere,

$$G^{sc}(r,r_0)=ik\sum_{lm}\left[s_M M_{lm}^{(h)}(r)\otimes M_{lm}^{(h,c)}(r_0)
+s_E N_{lm}^{(h)}(r)\otimes N_{lm}^{(h,c)}(r_0)\right].$$

The superscript c means angular conjugation only. It does not conjugate h_l(kr0). The dyad contracts the source index with p using an ordinary transpose. Physical work and intensity introduce p-dagger and E-star later.

At the local north-pole source, the electric radial m=0 harmonic is nonzero. The m=+/-1 surface gradients supply the two tangential directions; other m vanish at the source. Their exact pole limits are implemented to avoid dividing a vanishing associated Legendre function by sin(theta).

## 5. Total electromagnetic decay

The regular homogeneous Green tensor has Im Ghom=k I/(6pi). With vacuum classical peak-moment power

$$P_0=\frac{\omega^4p^\dagger p}{12\pi\epsilon_0c^3},$$

host-normalized decay is `FHost=1+(6pi/k) Im(p† Gsc p)/|p|²`. Performing the source-angle sums gives radial and tangential values

$$F_{tot,r}^{host}=1-\frac32\operatorname{Re}\sum_l(2l+1)l(l+1)a_l\left[\frac{h_l(x_0)}{x_0}\right]^2,$$
$$F_{tot,t}^{host}=1-\frac34\operatorname{Re}\sum_l(2l+1)
\left[a_l\left(\frac{\xi_l'(x_0)}{x_0}\right)^2+b_lh_l(x_0)^2\right].$$

No modulus square occurs in these feedback products. The Green wave leaves the source, scatters, and returns; its complex phase is crucial. The radial electric term contains a_l, while tangential feedback contains both electric and magnetic channels.

For arbitrary complex p define wr=|p_local,z|²/|p|². Spherical symmetry makes the coincident Green tensor diagonal in radial/tangential directions, so

$$F_{tot}^{host}=w_rF_{tot,r}^{host}+(1-w_r)F_{tot,t}^{host}.$$

The vacuum-normalized factor is nh times the host-normalized factor. Its quantum weak-coupling interpretation uses a transition moment and Gamma0=omega³|d|²/(3pi epsilon0 hbar c³); it is not a strong-coupling time-evolution model.

## 6. Radiative series and absorption

The outgoing field contains the original dipole and the induced sphere multipoles **coherently**. Angular orthogonality gives

$$F_{rad,r}^{host}=\frac32\sum_l(2l+1)l(l+1)
\left|\frac{j_l(x_0)-a_lh_l(x_0)}{x_0}\right|^2,$$
$$F_{rad,t}^{host}=\frac34\sum_l(2l+1)\left[
|j_l(x_0)-b_lh_l(x_0)|^2+
\left|\frac{\psi_l'(x_0)-a_l\xi_l'(x_0)}{x_0}\right|^2\right].$$

Weight by wr and multiply by nh for vacuum normalization. Even if the sphere has little contrast, enough multipoles must be retained to reconstruct the direct j_l source contribution: choosing order only from sphere size can underestimate the radiative series when kr0 is large.

For this **finite sphere** in a transparent host, no infinite planar guided channel exists, so Pabs=Ptotal-Prad (within numerical error). Lossless tests must give equality; passive absorbing spheres must give nonnegative absorption. Do not clamp a negative residual to zero to conceal a convergence problem.

## 7. Absolute far-field pattern

At large r, h_l(kr) approaches `(-i)^(l+1) exp(ikr)/(kr)`. Define source contractions

$$c^M_{lm}=s_Mh_l(x_0)\Phi_{lm}^*(\hat r_0)\cdot p,$$
$$c^N_{lm}=-s_Eh_l(x_0)\left[
H_0\Psi_{lm}^*(\hat r_0)\cdot p+
\frac{\sqrt{l(l+1)}}{x_0}Y_{lm}^*(\hat r_0)p_r\right],\quad H_0=\xi_l'(x_0)/\xi_l(x_0).$$

Using the dyadic expansion and N's tangential far limit, the scattered amplitude (coefficient of exp(ikr)/r) is

$$A_{sc}(\hat n)=\mu_0\omega^2 i\sum_{lm}(-i)^{l+1}
\left[\Phi_{lm}(\hat n)c^M_{lm}-i\Psi_{lm}(\hat n)c^N_{lm}\right].$$

Add the direct dipole amplitude

$$A_0(\hat n)=\frac{\mu_0\omega^2}{4\pi}
(I-\hat n\hat n^T)p\,e^{-ik\hat n\cdot r_0}.$$

Then

$$E_\infty(\hat n,r)=[A_0(\hat n)+A_{sc}(\hat n)]e^{ikr}/r,$$
$$\boxed{\frac{dP}{d\Omega}=\frac{n_h}{2\mu_0c}|A_0+A_{sc}|^2}.$$

All vectors are first computed in the local source-aligned frame and then rotated to global coordinates. Complex moment components interfere in a directional pattern even though their cross terms disappear from the integrated radial/tangential rates. Squaring the sphere field alone omits both direct emission and interference and cannot reproduce the decay formulas.

`mlsphere.radiative` integrates this absolute pattern over the full sphere. Comparing that result with the separate radiative multipole series tests the new angular implementation and its prefactor, not merely the normalized picture.

## 8. Stable modal products actually used by the code

Direct a_l times h_l(x0)^2 can underflow/overflow separately while their product remains finite. Define Hx=xi'/xi, Jx=psi'/psi and use the Wronskian

$$j_l(x)h_l(x)=\frac{i}{x^2(H_x-J_x)}.$$

Let rho_l=h_l(x0)/h_l(x). Then

$$s_{E/M}h_l(x_0)=K_{E/M}j_l(x)\rho_l,$$
$$s_{E/M}h_l(x_0)^2=K_{E/M}\,[j_l(x)h_l(x)]\,\rho_l^2.$$

The second expression is finite even when individual Hankel values are enormous. The code obtains J by downward logarithmic recurrence and H through outgoing Hankel ratios. For Q_l(z)=h_l(z)/h_(l-1)(z),

$$Q_1=1/z-i,\quad Q_l=(2l-1)/z-1/Q_{l-1},\quad H_l=1/Q_l-l/z,$$

and

$$\rho_l=\frac{x}{x_0}e^{i(x_0-x)}\prod_{j=1}^l Q_j(x_0)/Q_j(x).$$

The far-field coefficient needs s h, while total feedback needs s h²; both are stored directly. Standard a/b for the public Mie API are recovered by `a=-KE*j(x)^2/[j(x)h(x)]`, avoiding a separate huge h evaluation.

This is a double-precision normalized-product implementation derived from boundary matching and the Wronskian. It does not claim to reproduce every algorithm in the stable-Mie literature. Nonfinite results are rejected. Truncation is still a physical/numerical issue: increase nmax, check total/radiative balance and refine angular sampling separately. The hard input cap of 500 is not a certified accuracy range; extremely close emitters may require more capable specialized methods.

## 9. References and author code

1. C. F. Bohren and D. R. Huffman, *Absorption and Scattering of Light by Small Particles*, Wiley, [DOI](https://doi.org/10.1002/9783527618156). Conventional electric/magnetic Mie coefficients.
2. H. Mertens, A. F. Koenderink and A. Polman, “Plasmon-enhanced luminescence near noble-metal nanospheres: Comparison of exact theory and an improved Gersten and Nitzan model,” *Phys. Rev. B* **76**, 115123 (2007), [publication](https://doi.org/10.1103/PhysRevB.76.115123). This is the reference cited by the author's original project; the current derivation fixes its own coefficient convention explicitly.
3. M. Majic and E. C. Le Ru, “Numerically stable formulation of Mie theory for an emitter close to a sphere,” *Applied Optics* **59**, 1293–1300 (2020), [DOI](https://doi.org/10.1364/AO.379694). Motivation and broader numerical context for high-order overflow/underflow; this draft's measured range is reported in its own benchmarks.
4. [Author's original repository](https://github.com/knifelees3/DipoleEmissionNearNanoSphere). Executable original total-decay functions are preserved for direct comparison; ambiguous historical prose/figure labels are not used to redefine radial/tangential orientation.
