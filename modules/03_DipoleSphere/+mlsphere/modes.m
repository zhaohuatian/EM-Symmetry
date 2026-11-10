function B=modes(S)
%MODES Stable Mie products at exterior source, avoids tiny a_n * huge h_n^2.
% sE=-standard electric Mie a_n; sM=-standard magnetic Mie b_n.
N=S.nmax; l=1:N; n=S.nHost; ns=sqrt(S.epsilonSphere);
J=logj(S.x,N); Jy=logj(S.y,N); J0=logj(S.x0,N);
[H,Q]=logh(S.x,N); [H0,Q0]=logh(S.x0,N);
ratio=S.x/S.x0*exp(1i*(S.x0-S.x))*cumprod(Q0./Q);
prodjh=1i./(S.x^2*(H-J)); % Wronskian j(x)h(x), no large Hankel values
KE=(ns*J-n*Jy)./(n*Jy-ns*H);
KM=(ns*Jy-n*J)./(n*H-ns*Jy);
if S.epsilonSphere==S.nHost^2, KE(:)=0; KM(:)=0; end
ja=sqrt(pi/(2*S.x))*besselj(l+.5,S.x);
j0=sqrt(pi/(2*S.x0))*besselj(l+.5,S.x0);
B.SH_E=KE.*ja.*ratio; B.SH_M=KM.*ja.*ratio; % s_n h_n(x0)
B.SH2_E=KE.*prodjh.*ratio.^2; B.SH2_M=KM.*prodjh.*ratio.^2;
B.j=j0; B.jD=j0.*J0; B.H0=H0; B.l=l; B.S=S;
B.a=-KE.*ja.^2./prodjh; B.b=-KM.*ja.^2./prodjh;
assert(all(isfinite([B.SH_E B.SH_M B.SH2_E B.SH2_M B.jD])), ...
 'Numerical range exceeded; reduce nmax/size or use a specialized scaled solver.');
end
function L=logj(z,N)
% Downward recurrence of [z j_l(z)]'/[z j_l(z)].
K=N+ceil(abs(z))+60; a=0; L=zeros(1,N);
for j=K:-1:1
 a=j/z-1/(a+j/z);
 if j-1>=1 && j-1<=N, L(j-1)=a; end
end
end
function [L,Q]=logh(z,N)
Q=zeros(1,N); Q(1)=1/z-1i;
for j=2:N, Q(j)=(2*j-1)/z-1/Q(j-1); end
L=1./Q-(1:N)/z;
end
