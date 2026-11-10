function R=radiative(p,C,N)
%RADIATIVE Full-sphere power from absolute far-field integration.
if nargin<3, N=80; end
assert(N>=8&&N==floor(N),'N must be integer >=8.');
[theta,phi]=meshgrid(acos(-1+((1:N)-.5)*2/N),(0:2*N-1)*pi/N);
R.Far=mlsphere.farfield(p,C,theta,phi);
S=mlsphere.model(C); R.Prad=sum(R.Far.dPdOmega(:))*2*pi/N^2;
R.P0=S.omega^4*sum(abs(p).^2)/(12*pi*S.eps0*S.c^3);
R.Frad=R.Prad/R.P0;
end
