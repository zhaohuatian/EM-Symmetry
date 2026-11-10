# Derivation: interfaces, films and multilayers

## 1. Geometry, phasors and modal basis

Let layers be ordered along increasing z. The first and last are semi-infinite. For relative permittivity epsilon_l and nonmagnetic media,

$$k_0=2\pi/\lambda_0,\quad n_l=\sqrt{\epsilon_l},\quad k_l=k_0n_l,\quad q_l=\sqrt{k_l^2-\kappa^2}.$$

The time convention is exp(-i omega t). Select outgoing/decaying q: positive real for forward propagating waves, positive imaginary for evanescent waves. A common transverse wavevector kappa is conserved at all interfaces; do not reuse the same angle in materials of different index.

In the xz incidence plane define

$$s=(0,1,0)^T,\qquad p_l^+=(q_l,0,-\kappa)^T/k_l,\qquad p_l^-=(-q_l,0,-\kappa)^T/k_l.$$

An electric field in a layer is `s(A exp(iqz)+B exp(-iqz))` for TE and `p+ A exp(iqz)+p- B exp(-iqz)` for TM. Thus **TM Ex is proportional to A-B**, whereas its longitudinal Ez is proportional to -(A+B). This distinction is essential when comparing AB coefficient conventions.

Maxwell's equation gives H=(k cross E)/(mu0 omega) for each plane wave. Hence TE has Ey=A+B and Hx=-q(A-B)/(mu0 omega). TM has Ex=(q/k)(A-B) and Hy=n(A+B)/(mu0 c), after referencing both amplitudes to the same plane.

## 2. Single interface from boundary continuity

Set incident electric amplitude to 1, reflected to r and transmitted to t. No free surface current is applied. Tangential E/H must be continuous.

For TE:

$$1+r_s=t_s,\qquad q_i(1-r_s)=q_jt_s,$$

so

$$r_{ij}^s=\frac{q_i-q_j}{q_i+q_j},\qquad t_{ij}^s=\frac{2q_i}{q_i+q_j}.$$

For TM:

$$\frac{q_i}{k_i}(1-r_p)=\frac{q_j}{k_j}t_p,\qquad n_i(1+r_p)=n_jt_p,$$

which yields

$$r_{ij}^p=\frac{\epsilon_jq_i-\epsilon_iq_j}{\epsilon_jq_i+\epsilon_iq_j},\qquad t_{ij}^p=\frac{2n_in_jq_i}{\epsilon_jq_i+\epsilon_iq_j}.$$

These are unit **electric-polarization** amplitudes, not power-normalized or magnetic amplitudes. In this basis r_ji=-r_ij and t_ij t_ji=1-r_ij^2. The fact that the TM reflected basis changes its tangential sign explains why blindly comparing r_s and r_p at normal incidence can look inconsistent while the physical Cartesian reflected field is consistent.

## 3. Finite film: self-consistency and path sum

For media 1/2/3, film thickness d and f=exp(i q2 d), define A and B inside the film at its left boundary. Right incidence is zero. The equations are

$$B=r_{23}f^2A,\qquad A=t_{12}a+r_{21}B,$$

$$b=r_{12}a+t_{21}B,\qquad c=t_{23}fA.$$

Eliminating A/B gives

$$t_{\rm film}=\frac{t_{12}t_{23}f}{1-r_{21}r_{23}f^2},\quad
r_{\rm film}=r_{12}+\frac{t_{12}t_{21}r_{23}f^2}{1-r_{21}r_{23}f^2}.$$

Alternatively, each additional complete film round trip multiplies a path by g=r21 r23 f^2. Transmitted paths are `t12 t23 f (1+g+g^2+...)`; reflected paths are the direct r12 plus `t12 t21 r23 f^2 (1+g+...)`. For |g|<1 their sum equals the self-consistent solution. More generally the frequency-domain solution is defined by its outgoing/causal prescription; a perfectly closed resonant lossless cavity cannot sustain a finite driven steady state.

Substituting r21=-r12 recovers the familiar denominator 1+r12 r23 f^2. This sign change is bookkeeping, not a change of physics.

## 4. Two-port S matrix and cascade

Define incoming amplitudes [aL,aR] and outgoing [bL,bR]:

$$\begin{pmatrix}b_L\\b_R\end{pmatrix}=S\begin{pmatrix}a_L\\a_R\end{pmatrix}.$$

For an interface and a propagation segment,

$$S_{ij}=\begin{pmatrix}r_{ij}&t_{ji}\\t_{ij}&r_{ji}\end{pmatrix},\qquad
S_{\rm prop}=\begin{pmatrix}0&f\\f&0\end{pmatrix}.$$

The identity S is [0,1;1,0], not eye(2). For A followed by B, introduce connecting amplitudes u and v. They obey u=A21 aL+A22 v and v=B11 u+B12 aR. Eliminate them with Delta=1-A22 B11:

$$S_{11}=A_{11}+A_{12}B_{11}A_{21}/\Delta,\quad S_{12}=A_{12}B_{12}/\Delta,$$
$$S_{21}=B_{21}A_{21}/\Delta,\quad S_{22}=B_{22}+B_{21}A_{22}B_{12}/\Delta.$$

This scalar Redheffer product is `mlpw.cascade`. Each film inserts one propagation segment between two interfaces. It stores only decaying propagation factors in passive layers.

## 5. Transfer matrix and its range

From bR=S21 aL+S22 aR and bL=S11 aL+S12 aR,

$$\begin{pmatrix}a_L\\b_L\end{pmatrix}=
\begin{pmatrix}1/S_{21}&-S_{22}/S_{21}\\S_{11}/S_{21}&S_{12}-S_{11}S_{22}/S_{21}\end{pmatrix}
\begin{pmatrix}b_R\\a_R\end{pmatrix}.$$

Multiply these matrices in geometric order to obtain M. With unit left incidence, r=M21/M11 and t=1/M11. Its full S matrix has entries `[M21/M11, det(M)/M11; 1/M11, -M12/M11]`.

In an evanescent segment the transfer matrix contains 1/f, which grows exponentially. Double precision can overflow even though physical scattering remains finite. The code does not clip Im(qd), because clipping changes the physical problem. Use SMM when TMM is ill conditioned.

## 6. Stable internal fields

At each finite-layer left boundary, let P be the S matrix from the entrance to that plane, and Q the S matrix from that plane to the exit. For unit left incidence and zero right incidence,

$$A_{\rm left}=\frac{P_{21}}{1-P_{22}Q_{11}}.$$

At the right boundary, Qright is the remaining stack, so `B_right=Qright11 A_left exp(i q d)`. At distance x from the left boundary use

$$A(x)=A_{\rm left}e^{iqx},\qquad B(x)=B_{\rm right}e^{iq(d-x)}.$$

Both propagation distances are nonnegative. This avoids reconstructing an exponentially tiny B_left times an exponentially huge exp(-iqx), a combination that may become 0 times Inf numerically.

Construct E from s/p vectors and H from Maxwell's equation. The normal time-averaged flux is

$$S_z=\tfrac12\operatorname{Re}(E_xH_y^*-E_yH_x^*).$$

## 7. Amplitudes versus powers

For transparent propagating ports in this nonmagnetic unit-electric basis,

$$\mathcal R=|r|^2,\qquad \mathcal T=(q_{\rm out}/q_{\rm in})|t|^2,\qquad \mathcal A=1-\mathcal R-\mathcal T.$$

An evanescent port has no ordinary independent outgoing normal power in this formula. Outputs are NaN rather than misleading power fractions when its assumptions do not hold. A lossless stack must satisfy R+T=1. A reflection minimum alone does not establish perfect absorption.

## References and provenance

The author's [TransferMatrix source](https://github.com/knifelees3/TransferMatrix) is the implementation starting point; the complete two-port algebra is derived above. For thin-film internal-field applications, see C. A. Mack, “Analytical expression for the standing wave intensity in photoresist,” *Applied Optics* **25**, 1958–1961 (1986), [DOI](https://doi.org/10.1364/AO.25.001958). This draft's measured tests are described separately; it does not claim a new RCWA or COMSOL run.
