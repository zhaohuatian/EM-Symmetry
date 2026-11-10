function F=farfield(p,C,theta,phi,r)
%FARFIELD Total coherent direct+scattered exterior far field, SI units.
% theta/phi are GLOBAL spherical observation angles; 0<=theta<=pi.
if nargin<5, r=1; end
S=mlsphere.model(C); B=mlsphere.modes(S); p=p(:);
assert(numel(p)==3&&all(isfinite(p))&&norm(p)>0,'Invalid p.');
assert(isequal(size(theta),size(phi))&&all(theta(:)>=0&theta(:)<=pi),'Invalid observation angles.');
assert(isscalar(r)&&isfinite(r)&&r>0,'Invalid radius.');
ng=[sin(theta(:)).*cos(phi(:)),sin(theta(:)).*sin(phi(:)),cos(theta(:))];
n=ng*S.frame; t=acos(max(-1,min(1,n(:,3)))); f=atan2(n(:,2),n(:,1)); pl=S.frame.'*p;
direct=(pl.'-(n*pl).*n).*exp(-1i*S.x0*n(:,3))/(4*pi);
sc=zeros(size(direct));
for l=1:S.nmax
 for m=[0 1 -1]
  [Y,ps,ph]=mlsphere.angular(l,m,t,f);
  [Ys,pss,phs]=mlsphere.angular(l,m,0,0);
  cM=B.SH_M(l)*(conj(phs)*pl);
  cN=-B.SH_E(l)*(B.H0(l)*(conj(pss)*pl)+sqrt(l*(l+1))/S.x0*conj(Ys)*pl(3));
  sc=sc+1i*(-1i)^(l+1)*(ph*cM-1i*ps*cN);
 end
end
amplitude=S.mu0*S.omega^2*(direct+sc)*S.frame.';
assert(all(isfinite(amplitude(:))),'Angular series exceeded numerical range; inspect nmax and geometry.');
F.E=reshape(amplitude*exp(1i*S.k*r)/r,[size(theta) 3]);
F.dPdOmega=reshape(S.nHost/(2*S.mu0*S.c)*sum(abs(amplitude).^2,2),size(theta));
F.pattern=F.dPdOmega/max(F.dPdOmega(:)); F.theta=theta; F.phi=phi; F.r=r; F.nmax=S.nmax;
end
