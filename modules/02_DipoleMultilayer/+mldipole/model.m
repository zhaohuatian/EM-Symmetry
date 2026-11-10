function S=model(D)
% Validate original DinLayer geometry. SI units, bottom-to-top layer order.
S=D; S.ep=D.Eplist(:).'; S.z=D.dl(:).'; S.j=D.nDn+1;
assert(numel(S.ep)==numel(S.z)+1 && all(diff(S.z)>0),'Invalid interfaces.');
assert(S.j>=1 && S.j<=numel(S.ep),'Invalid source layer.');
assert(all(isfinite(S.ep)) && all(imag(S.ep)>=0),'Passive finite permittivity required.');
assert(imag(S.ep(S.j))==0 && real(S.ep(S.j))>0,'Source host must be lossless dielectric.');
assert(S.j==1 || D.POSD>S.z(S.j-1),'Source must be strictly inside its layer.');
assert(S.j==numel(S.ep) || D.POSD<S.z(S.j),'Source must be strictly inside its layer.');
assert(isfinite(D.WL0) && D.WL0>0 && isfinite(D.POSD),'Invalid wavelength/position.');
S.c=299792458; S.mu0=1.25663706212e-6; S.eps0=1/(S.mu0*S.c^2);
S.k0=2*pi/D.WL0; S.omega=S.c*S.k0; S.k=S.k0*sqrt(S.ep);
S.nh=sqrt(S.ep(S.j));
end
