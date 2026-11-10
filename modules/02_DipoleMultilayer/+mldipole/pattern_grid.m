function [pat,F]=pattern_grid(p,D,objective,sgn)
% Original kx_grid/ky_grid interface. Spherical W/sr (objective=0), or
% dP/d(ux)d(uy) in W for ux=kx/k_out, uy=ky/k_out (objective=1).
% Outside the exterior light cone: NaN. No calibrated objective model.
if nargin<3, objective=0; end
if nargin<4, sgn=1; end
assert(ismember(objective,[0 1]),'objective must be 0 or 1.');
S=mldipole.model(D); o=1+(sgn>0)*(numel(S.ep)-1);
assert(isreal(S.k(o)) && S.k(o)>0,'Lossless exterior required.');
rho=hypot(D.kx_grid,D.ky_grid); mask=rho<S.k(o);
theta=asin(rho(mask)/S.k(o)); phi=atan2(D.ky_grid(mask),D.kx_grid(mask));
F=mldipole.farfield(p,D,theta,phi,sgn);
pat=nan(size(rho)); pat(mask)=F.dPdOmega./cos(theta).^objective;
end
