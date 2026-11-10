function S=diagonal_stack(C,kx,ky,Ain,Bin)
%DIAGONAL_STACK Four-mode boundary solve for a reciprocal diagonal tensor stack.
% C.lambda [m], C.epsilonDiagonal Lx3, C.interfaces increasing [m].
% Ain/Bin are 2xK incident TE/TM amplitudes at bottom/top interface planes.
% Every finite layer has two upward and two downward eigenwaves. Unlike a
% transfer product, opposite-going waves are anchored to opposite boundaries.
if nargin<4, Ain=eye(2); Bin=zeros(2); end
assert(isscalar(C.lambda)&&isreal(C.lambda)&&isfinite(C.lambda)&&C.lambda>0,'Positive finite vacuum wavelength required.');
ep=C.epsilonDiagonal; z=C.interfaces(:).'; L=size(ep,1); k0=2*pi/C.lambda;
assert(size(ep,2)==3&&L>=2&&numel(z)==L-1&&all(diff(z)>0),'Invalid tensor geometry.');
assert(all(isfinite(ep(:)))&&all(imag(ep(:))>=0)&&all(abs(ep(:))>0),'Passive diagonal nonzero epsilon required.');
assert(isequal(size(Ain),size(Bin))&&size(Ain,1)==2,'Incoming amplitude matrices must match.');
assert(all(ep(1,:)==ep(1,1))&&all(ep(end,:)==ep(end,1)),'Exterior layers must be isotropic.');
u=kx/k0; v=ky/k0; modes=cell(1,L); qs=cell(1,L); lo=zeros(1,L); hi=lo;
for l=1:L
 if l==1, lo(l)=z(1); hi(l)=z(1);
 elseif l==L, lo(l)=z(end); hi(l)=z(end);
 else, lo(l)=z(l-1); hi(l)=z(l); end
 if all(ep(l,:)==ep(l,1))
  n=sqrt(ep(l,1)); w=sqrt(ep(l,1)-u^2-v^2); if imag(w)<0, w=-w; end
  assert(abs(w)>1e-11,'Avoid an exact modal light line.');
  ph=atan2(v,u); sv=[-sin(ph);cos(ph);0]; rho=sqrt(u^2+v^2);
  E=[sv,[w*cos(ph);w*sin(ph);-rho]/n,sv,[-w*cos(ph);-w*sin(ph);-rho]/n];
  K=[u u u u;v v v v;w w -w -w]; H=cross(K,E,1);
  modes{l}=[E(1:2,:);H(1:2,:)]; qs{l}=[w w -w -w];
 else
  ex=ep(l,1); ey=ep(l,2); ez=ep(l,3);
  M=[0 0 u*v/ez 1-u^2/ez;0 0 v^2/ez-1 -u*v/ez; ...
      -u*v u^2-ey 0 0;ex-v^2 u*v 0 0];
  [V,Q]=eig(M); q=diag(Q).'; up=false(1,4);
  for j=1:4
   if abs(imag(q(j)))>1e-9, up(j)=imag(q(j))>0;
   else, flux=real(V(1,j)*conj(V(4,j))-V(2,j)*conj(V(3,j))); up(j)=flux>0; end
  end
  assert(sum(up)==2,'Cannot separate two outgoing/incoming modes for this tensor.');
  ids=[find(up),find(~up)]; modes{l}=V(:,ids); qs{l}=q(ids);
 end
end
% Tangential state continuity at every physical interface, and exterior inputs.
M=zeros(4*L); rhs=zeros(4*L,size(Ain,2));
M(1:2,1:2)=eye(2); rhs(1:2,:)=Ain;
M(3:4,4*L-1:4*L)=eye(2); rhs(3:4,:)=Bin;
for l=1:L-1
 rows=4*l+(1:4);
 M(rows,4*l-3:4*l)=at(l,z(l));
 M(rows,4*l+1:4*l+4)=-at(l+1,z(l));
end
coeff=M\rhs;
S.C=C; S.k0=k0; S.kx=kx; S.ky=ky; S.modes=modes; S.q=qs;
S.lowerReference=lo; S.upperReference=hi; S.coefficients=reshape(coeff,[4,L,size(rhs,2)]);
S.residual=norm(M*coeff-rhs,'fro')/max(1,norm(rhs,'fro'));
assert(isfinite(S.residual)&&S.residual<1e-7,'Boundary solve did not meet residual tolerance.');
    function V=at(l,zz)
        phases=[exp(1i*k0*qs{l}(1:2)*(zz-lo(l))),exp(1i*k0*qs{l}(3:4)*(zz-hi(l)))];
        V=modes{l}.*phases;
    end
end
