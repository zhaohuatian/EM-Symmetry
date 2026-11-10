function S=model(C)
%MODEL Validate an EXTERIOR dipole and homogeneous nonmagnetic sphere.
% lambda/radius/position in m; epsilonSphere relative; nHost real >0.
assert(isscalar(C.lambda)&&isfinite(C.lambda)&&C.lambda>0,'Invalid lambda.');
assert(isscalar(C.radius)&&isfinite(C.radius)&&C.radius>0,'Invalid radius.');
assert(isscalar(C.nHost)&&isreal(C.nHost)&&C.nHost>0,'Lossless positive host required.');
assert(isscalar(C.epsilonSphere)&&isfinite(C.epsilonSphere)&&abs(C.epsilonSphere)>0 ...
 && imag(C.epsilonSphere)>=0,'Passive nonzero sphere epsilon required.');
S=C; S.position=C.position(:);
assert(numel(S.position)==3&&isreal(S.position)&&all(isfinite(S.position)),'Invalid source position.');
S.r0=norm(S.position); assert(S.r0>C.radius,'Source must be outside the sphere.');
S.c=299792458; S.mu0=1.25663706212e-6; S.eps0=1/(S.mu0*S.c^2);
S.k0=2*pi/C.lambda; S.k=S.k0*C.nHost; S.omega=S.k0*S.c;
S.x=S.k*C.radius; S.x0=S.k*S.r0; S.y=S.k0*sqrt(C.epsilonSphere)*C.radius;
if ~isfield(C,'nmax')
 S.nmax=max([ceil(S.x+4*S.x^(1/3)+8),ceil(S.x0+4*S.x0^(1/3)+8), ...
     ceil(log(1e-10)/(2*log(C.radius/S.r0)))]);
end
assert(isscalar(S.nmax)&&S.nmax>=1&&S.nmax==floor(S.nmax)&&S.nmax<=500,'nmax must be 1..500; check convergence.');
er=S.position/S.r0; [~,j]=min(abs(er)); a=zeros(3,1); a(j)=1;
ex=cross(a,er); ex=ex/norm(ex); ey=cross(er,ex); S.frame=[ex ey er];
end
