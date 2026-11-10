function F=modal_field(S,z)
%MODAL_FIELD All-layer plane-wave E/H from diagonal_stack's A/B coefficients.
% No dipole source in this auxiliary solution. F.E is N x 3 x K excitations.
z=z(:); K=size(S.coefficients,3); E=zeros(numel(z),3,K); H=E;
mu0=1.25663706212e-6; c=299792458; u=S.kx/S.k0; v=S.ky/S.k0;
for a=1:numel(z)
 l=1+sum(z(a)>=S.C.interfaces(:));
 q=S.q{l}; phase=[exp(1i*S.k0*q(1:2)*(z(a)-S.lowerReference(l))), ...
    exp(1i*S.k0*q(3:4)*(z(a)-S.upperReference(l)))];
 state=(S.modes{l}.*phase)*reshape(S.coefficients(:,l,:),4,K);
 Ez=(v*state(3,:)-u*state(4,:))/S.C.epsilonDiagonal(l,3);
 hz=u*state(2,:)-v*state(1,:);
 E(a,:,:)=reshape([state(1:2,:);Ez],1,3,K);
 H(a,:,:)=reshape([state(3:4,:);hz]/(mu0*c),1,3,K);
end
F.z=z; F.E=E; F.H=H;
end
