function F=field(R,z)
%FIELD Cartesian E [incident amplitude units], H [A/m for unit V/m].
% Layer coordinate increases from left to right; incidence plane xz.
% At an interface the sample belongs to the layer on its right.
z=z(:); e=R.C.epsilon(:).'; L=numel(e); c=299792458; mu0=1.25663706212e-6;
E=zeros(numel(z),3); H=E; layer=zeros(size(z));
for a=1:numel(z)
 l=1+sum(z(a)>=R.z); layer(a)=l;
 if l==1
  A=exp(1i*R.q(l)*z(a)); B=R.Bright(l)*exp(-1i*R.q(l)*z(a));
 elseif l==L
  A=R.Aleft(l)*exp(1i*R.q(l)*(z(a)-R.z(end))); B=0;
 else
  x=z(a)-R.z(l-1); d=R.C.thickness(l-1);
  A=R.Aleft(l)*exp(1i*R.q(l)*x);
  B=R.Bright(l)*exp(1i*R.q(l)*(d-x));
 end
 n=sqrt(e(l)); co=R.q(l)/(R.k0*n); si=R.C.kparallel/(R.k0*n);
 if strcmpi(R.C.polarization,'TE')
  E(a,:)=[0 A+B 0]; H(a,:)=n/(mu0*c)*[-co*(A-B) 0 si*(A+B)];
 else
  E(a,:)=[co*(A-B) 0 -si*(A+B)]; H(a,:)=n/(mu0*c)*[0 A+B 0];
 end
end
F.z=z; F.E=E; F.H=H; F.layer=layer;
F.Sz=.5*real(F.E(:,1).*conj(F.H(:,2))-F.E(:,2).*conj(F.H(:,1)));
end
