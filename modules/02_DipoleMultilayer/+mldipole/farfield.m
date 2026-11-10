function F=farfield(p,D,theta,phi,sgn,r)
% p: complex peak electric moment [C m]. theta: angle from outward +/-z.
% F.E [V/m], F.dPdOmega [W/sr], F.pattern normalized only for display.
if nargin<5, sgn=1; end
if nargin<6, r=1; end
S=mldipole.model(D); p=p(:); assert(numel(p)==3 && all(isfinite(p)) && norm(p)>0,'Invalid p.');
assert(ismember(sgn,[-1 1]) && isscalar(r) && r>0,'Invalid side/radius.');
assert(isequal(size(theta),size(phi)) && all(theta(:)>=0 & theta(:)<pi/2),'Angles must match and 0<=theta<pi/2.');
o=1+(sgn>0)*(numel(S.ep)-1); no=sqrt(S.ep(o));
assert(isreal(no) && no>0,'Far-field exterior must be lossless dielectric.');
E=zeros(numel(theta),3); C=S.mu0*S.omega^2/(4*pi);
for a=1:numel(theta)
    rho=S.k(o)*sin(theta(a)); ph=phi(a);
    [Ru,Rd,Tu,Td,q]=mldipole.rt(S,rho);
    sv=[-sin(ph);cos(ph);0];
    pp=[q(S.j)*cos(ph);q(S.j)*sin(ph);-rho]/S.k(S.j);
    pm=[-q(S.j)*cos(ph);-q(S.j)*sin(ph);-rho]/S.k(S.j);
    ap=[sv.'*p,pp.'*p]; am=[sv.'*p,pm.'*p];
    den=1-Ru.*Rd;
    if sgn>0, b=Tu.*(ap+Rd.*am)./den; else, b=Td.*(am+Ru.*ap)./den; end
    po=[sgn*q(o)*cos(ph);sgn*q(o)*sin(ph);-rho]/S.k(o);
    if o==S.j, zb=S.POSD; elseif sgn>0, zb=S.z(end); else, zb=S.z(1); end
    v=C*q(o)/q(S.j)*(sv*b(1)+po*b(2))*exp(-1i*sgn*q(o)*zb);
    E(a,:)=v.'*exp(1i*S.k(o)*r)/r;
end
F.E=reshape(E,[size(theta),3]);
F.dPdOmega=reshape(no/(2*S.mu0*S.c)*sum(abs(E).^2,2)*r^2,size(theta));
F.pattern=F.dPdOmega/max(F.dPdOmega(:));
F.theta=theta; F.phi=phi; F.radius=r; F.side=sgn;
end
