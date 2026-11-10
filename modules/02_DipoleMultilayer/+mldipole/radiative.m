function R=radiative(p,D,N)
% Full two-exterior hemispheres. Midpoint in cos(theta), periodic azimuth.
if nargin<3, N=160; end
assert(N>=8 && N==floor(N),'N must be an integer >=8.');
mu=((1:N)-0.5)/N; phi=(0:2*N-1)*pi/N;
[tt,pp]=meshgrid(acos(mu),phi);
R.Up=mldipole.farfield(p,D,tt,pp,1);
R.Dn=mldipole.farfield(p,D,tt,pp,-1);
R.PUp=sum(R.Up.dPdOmega(:))*pi/N^2;
R.PDn=sum(R.Dn.dPdOmega(:))*pi/N^2;
S=mldipole.model(D);
R.P0=S.omega^4*sum(abs(p).^2)/(12*pi*S.eps0*S.c^3);
R.Prad=R.PUp+R.PDn; R.Frad=R.Prad/R.P0;
end
