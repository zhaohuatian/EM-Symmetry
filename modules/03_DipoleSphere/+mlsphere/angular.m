function [Y,Psi,Phi]=angular(l,m,theta,phi)
%ANGULAR Normalized Y_lm, surface gradient/sqrt(l(l+1)), rhat cross Psi.
% Only m=0,+/-1 needed with source rotated to north pole; exact pole limits.
t=theta(:); f=phi(:); ct=cos(t); st=sin(t); N=numel(t);
P=legendre(l,ct.'); P0=P(1,:).'; P1=P(2,:).';
if m==0
 v=sqrt((2*l+1)/(4*pi)); Y=v*P0; dt=v*P1; az=zeros(N,1);
else
 v=exp(.5*(log((2*l+1)/(4*pi))+gammaln(l)-gammaln(l+2)));
 if l==1, prev=zeros(N,1); else, a=legendre(l-1,ct.'); prev=a(2,:).'; end
 d=(l*ct.*P1-(l+1)*prev)./st; ratio=P1./st;
 north=abs(t)<1e-10; south=abs(t-pi)<1e-10;
 aa=l*(l+1)/2; d(north)=-aa; ratio(north)=-aa;
 d(south)=(-1)^(l+1)*aa; ratio(south)=(-1)^l*aa;
 Y=v*P1.*exp(1i*f); dt=v*d.*exp(1i*f); az=1i*v*ratio.*exp(1i*f);
 if m==-1, Y=-conj(Y); dt=-conj(dt); az=-conj(az); end
end
et=[ct.*cos(f),ct.*sin(f),-st]; ef=[-sin(f),cos(f),zeros(N,1)];
Psi=(dt.*et+az.*ef)/sqrt(l*(l+1));
Phi=(dt.*ef-az.*et)/sqrt(l*(l+1));
end
