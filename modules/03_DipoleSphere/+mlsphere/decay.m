function R=decay(p,C)
%DECAY Total, radiative and absorption enhancement for arbitrary complex p.
% F* normalized to vacuum; F*Host normalized to homogeneous host.
S=mlsphere.model(C); B=mlsphere.modes(S); p=p(:);
assert(numel(p)==3&&all(isfinite(p))&&norm(p)>0,'Invalid dipole moment.');
pl=S.frame.'*p; w=abs(pl(3))^2/sum(abs(pl).^2); l=B.l; x=S.x0;
tr=1+1.5*real(sum((2*l+1).*l.*(l+1).*B.SH2_E/x^2));
tt=1+.75*real(sum((2*l+1).*(B.SH2_E.*B.H0.^2+B.SH2_M)));
rr=1.5*sum((2*l+1).*l.*(l+1).*abs((B.j+B.SH_E)/x).^2);
rt=.75*sum((2*l+1).*(abs(B.j+B.SH_M).^2+abs(B.jD+B.SH_E.*B.H0).^2));
R.FtotalHost=w*tr+(1-w)*tt; R.FradHost=w*rr+(1-w)*rt;
R.Ftotal=S.nHost*R.FtotalHost; R.Frad=S.nHost*R.FradHost;
R.Fabs=R.Ftotal-R.Frad; R.radial=struct('totalHost',tr,'radHost',rr);
R.tangential=struct('totalHost',tt,'radHost',rt);
R.P0=S.omega^4*sum(abs(p).^2)/(12*pi*S.eps0*S.c^3);
R.Ptotal=R.P0*R.Ftotal; R.Prad=R.P0*R.Frad; R.Pabs=R.Ptotal-R.Prad;
R.nmax=S.nmax; R.Gamma0=S.omega^3*sum(abs(p).^2)/(3*pi*S.eps0*1.054571817e-34*S.c^3);
R.Gamma=R.Gamma0*R.Ftotal;
end
