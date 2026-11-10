function F=reciprocal_farfield(p,C,theta,phi,side,r)
%RECIPROCAL_FARFIELD Absolute far field for diagonal anisotropic layers.
% C: lambda, epsilonDiagonal Lx3, interfaces, sourceZ.
% Far-field ONLY: no anisotropic total self-decay implementation is implied.
if nargin<5, side=-1; end
if nargin<6, r=1; end
p=p(:); assert(numel(p)==3&&all(isfinite(p))&&norm(p)>0,'Invalid p.');
assert(isequal(size(theta),size(phi))&&all(theta(:)>=0&theta(:)<pi/2),'Invalid angles.');
assert(ismember(side,[-1 1])&&isscalar(r)&&r>0,'Invalid side/radius.');
L=size(C.epsilonDiagonal,1); o=1+(side>0)*(L-1); eo=C.epsilonDiagonal(o,1);
assert(isreal(eo)&&eo>0,'Transparent isotropic observation exterior required.');
assert(~any(C.sourceZ==C.interfaces),'Source must be strictly inside a layer.');
no=sqrt(eo); mu0=1.25663706212e-6; c=299792458; k0=2*pi/C.lambda;
E=zeros(numel(theta),3); residual=zeros(size(theta));
for j=1:numel(theta)
 t=theta(j); ph=phi(j); rho=k0*no*sin(t);
 kx=-rho*cos(ph); ky=-rho*sin(ph);
 % Desired reciprocal electric polarization equals the outgoing basis.
 % For rho=0 atan2(0,0) is basis-arbitrary: project explicitly instead of
 % relying on a fixed TE sign at normal incidence.
 es=[-sin(ph);cos(ph);0]; ep=[side*cos(t)*cos(ph);side*cos(t)*sin(ph);-sin(t)];
 pp=atan2(ky,kx); sr=[-sin(pp);cos(pp);0];
 pr=[-side*cos(t)*cos(pp);-side*cos(t)*sin(pp);-sin(t)];
 inc=[sr.'*es,sr.'*ep;pr.'*es,pr.'*ep];
 Ain=zeros(2); Bin=Ain; if side>0, Bin=inc; else, Ain=inc; end
 S=mlpw.diagonal_stack(C,kx,ky,Ain,Bin);
 A=mlpw.modal_field(S,C.sourceZ); aux=reshape(A.E,3,2);
 amp=mu0*(c*k0)^2/(4*pi)*(es*(aux(:,1).'*p)+ep*(aux(:,2).'*p));
 zb=C.interfaces(1); if side>0, zb=C.interfaces(end); end
 amp=amp*exp(-1i*side*k0*no*cos(t)*zb);
 E(j,:)=amp.'*exp(1i*k0*no*r)/r; residual(j)=S.residual;
end
F.E=reshape(E,[size(theta) 3]);
F.dPdOmega=reshape(no/(2*mu0*c)*sum(abs(E).^2,2)*r^2,size(theta));
F.pattern=F.dPdOmega/max(F.dPdOmega(:)); F.theta=theta; F.phi=phi;
F.side=side; F.r=r; F.maxBoundaryResidual=max(residual(:));
end

