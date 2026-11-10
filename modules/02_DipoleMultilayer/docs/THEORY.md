# Derivation: angular-spectrum source, layered feedback and decay

For a guided reading route with every-layer A+B fields, interface matching, diagrams and examples, start with the [worked tutorial](TUTORIAL.md). The [diagonal-anisotropic extension](FIGURE2_RECONSTRUCTION.md) has a separate four-mode derivation; the scalar formulas below refer to isotropic layers.

## 1. Geometry and source normalization

Layers increase from bottom to top. A dipole p in source layer h is at r0=(0,0,z0). Let epsilon_l be relative permittivity and mu=mu0. Use peak phasors exp(-i omega t), k0=omega/c and k_l=k0 sqrt(epsilon_l).

The polarization source is P=p delta(r-r0), so J=-i omega p delta(r-r0). Maxwell's equations yield

$$[\nabla\times\nabla\times-k_0^2\epsilon(r)]G(r,r_0)=I\delta(r-r_0),\qquad
E=\mu_0\omega^2Gp.$$

The source factor is omega² for a dipole in C m. `i omega mu0 G` instead acts on a **current moment**. Mixing the two definitions preserves a normalized image but destroys absolute field and power values.

For a uniform source host,

$$G_h^0=(I+\nabla\nabla/k_h^2)\frac{e^{ik_h|r-r_0|}}{4\pi|r-r_0|}.$$

Its singular real self-part must not be evaluated as an ordinary finite field at r=r0.

## 2. Weyl representation with an explicit Fourier convention

Use `G=(2pi)^-2 integral Gtilde exp(i kparallel.rho) d²kparallel`. Write kparallel=kappa(cos phi,sin phi), q_l=sqrt(k_l²-kappa²), and select the outgoing/decaying branch. The scalar transformed Green function solves a one-dimensional derivative jump and is

$$\widetilde g_h=\frac{i}{2q_h}e^{iq_h|z-z_0|}.$$

Define TE and TM vectors

$$s=(-\sin\phi,\cos\phi,0)^T,\quad
e_l^\pm=(\pm q_l\cos\phi,\pm q_l\sin\phi,-\kappa)^T/k_l.$$

Then `s s^T + e± e±^T = I-k± k±^T/k²`. For complex wavevectors this is a bilinear transpose identity; do not conjugate source projections. At separated points the host dyadic spectrum is

$$\widetilde G_h^0=\frac{i}{2q_h}(ss^T+e_h^\pm(e_h^\pm)^T)e^{iq_h|z-z_0|}.$$

The full distribution also contains a longitudinal contact term, which is not needed in separated far fields or the finite scattered self-field.

## 3. Scatter each TE/TM angular component

For each kappa independently, use the interface coefficients derived in the plane-wave module:

$$r_s=(q_i-q_j)/(q_i+q_j),\quad
r_p=(\epsilon_jq_i-\epsilon_iq_j)/(\epsilon_jq_i+\epsilon_iq_j),$$
$$t_s=2q_i/(q_i+q_j),\quad t_p=2n_in_jq_i/(\epsilon_jq_i+\epsilon_iq_j).$$

Starting from the exterior R=0,T=1, add an interface and the neighboring finite layer with f=exp(i q_j d_j):

$$R_i=\frac{r_{ij}+R_jf^2}{1+r_{ij}R_jf^2},\qquad
T_i=\frac{t_{ij}T_jf}{1+r_{ij}R_jf^2}.$$

This is the first column of the two-port S cascade, not an approximation. The reflected coefficient is referenced to the current entrance plane. Finally move the reference plane from the source-layer interface to the dipole plane:

$$R_u=\bar R_ue^{2iq_hd_u},\quad T_u=\bar T_ue^{iq_hd_u},\quad
R_d=\bar R_de^{2iq_hd_d},\quad T_d=\bar T_de^{iq_hd_d}.$$

These Ru/Rd already contain the propagation to each stack and back. A missing side has R=0,T=1. All interface-internal reflection is included in these stack responses; the remaining round trips between the two stacks are solved next.

## 4. Dipole launches and two-sided feedback

After factoring out i mu0 omega²/(2q_h), the directly emitted amplitudes are

$$a_s^+=a_s^-=s^Tp,\qquad a_p^\pm=(e_h^\pm)^Tp.$$

Explicitly a_s=-px sin phi+py cos phi and a_p±=±(q_h/k_h)(px cos phi+py sin phi)-(kappa/k_h)pz. Horizontal and vertical source components therefore have different signs on TM reversal.

At the **source plane**, the total upward/downward amplitudes satisfy

$$U=a^++R_dV,\qquad V=a^-+R_uU.$$

The source plane is a virtual reference inside a homogeneous host, not another dielectric interface. A reflected return crosses it freely while the maintained source adds its direct emission. Upward incidence on the upper stack produces a returned downward wave Ru U **and** an external transmitted wave Tu U; only the former enters the feedback equation.

Substitution gives `(1-Rd Ru)U=a+ + Rd a-`. Thus

$$D=1-R_uR_d,\quad U=\frac{a^++R_da^-}{D},\quad V=\frac{a^-+R_ua^+}{D}.$$

D is simply the determinant of `[1,-Rd;-Ru,1]`. The path interpretation, when convergent, is

$$U=(1+g+g^2+\cdots)(a^++R_da^-),\quad
V=(1+g+g^2+\cdots)(a^-+R_ua^+),\quad g=R_uR_d.$$

Each extra g is a complete upper/lower round trip. It is 1/D, not D, that expands as this geometric series. A causal outgoing solution is needed at lossless modal poles; the ordinary series cannot be assumed convergent for every evanescent angular component.

## 5. Absolute exterior far field

For upper/lower exterior o, set kappa=k_o sin theta, q_o=k_o cos theta. The source denominator remains q_h. Let b_s/p+=Tu_s/p U_s/p and b_s/p-=Td_s/p V_s/p. If the last exterior interface is at zb, the spectral field in an inverse transform with no outside prefactor is

$$A^\pm=\frac{i\mu_0\omega^2}{8\pi^2q_h}
e^{\mp iq_oz_b}(s b_s^\pm+e_o^\pm b_p^\pm).$$

Two-dimensional stationary phase contributes `-2pi i q_o exp(i k_o r)/r`, giving

$$\boxed{E_\infty^\pm=\frac{\mu_0\omega^2}{4\pi r}e^{ik_or}
\frac{q_o}{q_h}e^{\mp iq_oz_b}(s b_s^\pm+e_o^\pm b_p^\pm)}.$$

If source and exterior coincide without an intervening interface, zb=z0. The reference-plane phase matters for complex-field comparisons. The formula is a 1/r asymptote, not a complete near-field solution. A source-layer-evanescent component can still radiate into a higher-index substrate, so use the **observation** light cone to set kappa.

For a homogeneous host this reduces to

$$E_\infty=\frac{\mu_0\omega^2}{4\pi r}e^{ikr-ik\hat n\cdot r_0}(I-\hat n\hat n^T)p,$$

which supplies an absolute amplitude and phase check independent of normalized patterns.

## 6. Angular power and radiative enhancement

In a lossless exterior, H=(n_o/mu0 c) n_hat cross E, so

$$I_o=\frac{dP_o}{d\Omega}=\frac{n_or^2}{2\mu_0c}|E_\infty|^2,\quad
P_{\rm rad}=\sum_{u,d}\int_0^{2\pi}\int_0^{\pi/2}I_o\sin\theta\,d\theta d\phi.$$

Vacuum reference power and radiative factor are

$$P_0=\frac{\omega^4p^\dagger p}{12\pi\epsilon_0c^3},\qquad F_{\rm rad}=P_{\rm rad}/P_0.$$

A maximum-normalized pattern has no absolute power information. For pupil coordinates ux=kx/kout,uy=ky/kout, `dP/dux duy=I/cos(theta)`. This is a geometric Jacobian, not a full calibrated objective/detector model.

## 7. Local scattered Green tensor

Write G=Ghom+Gsc. For a transparent source host, `Im Ghom(r0,r0)=k_h I/(6pi)`. The positive power provided by the maintaining source is

$$P_{\rm total}=\frac\omega2\operatorname{Im}(p^\dagger E_{\rm reg})
=n_hP_0+\frac{\mu_0\omega^3}{2}\operatorname{Im}(p^\dagger G^{sc}p).$$

For a planar isotropic stack the coincident scattered dyadic is diagonal with xx=yy=Gparallel and zz=Gperp. Define

$$A_\sigma=\frac{R_u^\sigma+R_d^\sigma+2R_u^\sigma R_d^\sigma}{1-R_u^\sigma R_d^\sigma},\quad
B_p=\frac{R_u^p+R_d^p-2R_u^pR_d^p}{1-R_u^pR_d^p}.$$

The angular average of reflected TE/TM dyads gives

$$G_\parallel^{sc}=\frac{i}{8\pi}\int_0^\infty\frac\kappa{q_h}
\left(A_s-\frac{q_h^2}{k_h^2}B_p\right)d\kappa,$$
$$G_\perp^{sc}=\frac{i}{4\pi}\int_0^\infty\frac{\kappa^3}{q_hk_h^2}A_p\,d\kappa.$$

The minus sign in the horizontal TM contribution follows from the opposite tangential directions of e+ and e-. Both vertical components have the same sign. With wperp=|pz|²/|p|² and wparallel=1-wperp,

$$F_{\rm total}=n_h+\frac{6\pi}{k_0}
\left[w_\parallel\operatorname{Im}G_\parallel^{sc}+w_\perp\operatorname{Im}G_\perp^{sc}\right].$$

Integrated decay depends only on radial weights relative to the plane, while directional patterns generally retain complex-polarization interference.

## 8. Numerical contour, normalization and limits

The kernel is integrated using u=kappa/k0=t-i eta t/(1+t), derivative `1-i eta/(1+t)^2`, followed by a straight endpoint segment back to the real cutoff. Both Jacobians are included. The default source-distance-based cutoff is a starting point: refine cutoff, eta and quadrature tolerances. No artificial material loss is introduced. Very narrow modes and extreme parameter ranges require dedicated checks of analytic continuation and convergence.

For infinite homogeneous index nh, Ftotal=Frad=nh, without microscopic local-field correction. Divide both factors by nh for host-normalized results. For a quantum transition dipole d, Gamma0=omega³|d|²/(3pi epsilon0 hbar c³) and Gamma/Gamma0=Ftotal. If one uses the same numerical value for a classical peak moment p=d, Gamma=4Ptotal/(hbar omega); the factor four is an amplitude convention, not an enhancement correction.

`Ftotal-Frad` includes guided and absorbed channels in a planar stack. It is not solely absorption. Lossy source hosts, source-on-interface limits and strong-coupling dynamics are outside this implementation.

## References

1. M. Paulus, P. Gay-Balmaz, O. J. F. Martin, “Accurate and efficient computation of the Green's tensor for stratified media,” *Phys. Rev. E* **62**, 5797–5807 (2000), [DOI](https://doi.org/10.1103/PhysRevE.62.5797).
2. J. E. Sipe, “New Green-function formalism for surface optics,” *JOSA B* **4**, 481–489 (1987), [DOI](https://doi.org/10.1364/JOSAB.4.000481).
3. C. Creatore and L. C. Andreani, “Quantum theory of spontaneous emission in multilayer dielectric structures,” *Phys. Rev. A* **78**, 063825 (2008), [open manuscript](https://arxiv.org/abs/0810.4473).
4. Original author code and notes are identified in [provenance](../PROVENANCE.md). Equations here explicitly fix their own source/phase conventions; do not mix unexplained prefactors from another Green definition.
