function R=solve(C,method)
%SOLVE Complex reflection/transmission and stable internal coefficients.
% C.lambda [m], C.epsilon [left,...,right], C.thickness [m] for finite
% layers only, C.kparallel [1/m], C.polarization ('TE'/'TM').
% method: 'smm' (default) or 'tmm' (no phase clipping).
if nargin<2, method='smm'; end
assert(ismember(lower(string(method)),["smm","tmm"]),'Unknown method.');
e=C.epsilon(:).'; d=C.thickness(:).'; L=numel(e);
assert(L>=2 && numel(d)==L-2 && all(isfinite(d)) && all(d>0),'Invalid layers/thickness.');
assert(all(isfinite(e)) && all(imag(e)>=0) && all(abs(e)>0),'Passive nonzero epsilon required.');
assert(isscalar(C.lambda) && C.lambda>0 && isfinite(C.lambda),'Invalid wavelength.');
assert(isscalar(C.kparallel) && isfinite(C.kparallel),'Invalid kparallel.');
k0=2*pi/C.lambda; q=sqrt(k0^2*e-C.kparallel^2); q(imag(q)<0)=-q(imag(q)<0);
assert(all(abs(q)>1e-13*k0),'Avoid exact light-line points; take a limit.');
% Segment chain and reference-plane locations for each finite layer.
seg={}; left=zeros(1,L); right=left;
for j=1:L-1
 seg{end+1}=mlpw.interface(e(j),e(j+1),q(j),q(j+1),C.polarization); %#ok<AGROW>
 left(j+1)=numel(seg);
 if j+1<L
  f=exp(1i*q(j+1)*d(j)); seg{end+1}=[0 f;f 0]; %#ok<AGROW>
 end
 right(j+1)=numel(seg);
end
J=[0 1;1 0]; M=numel(seg); pre=cell(1,M+1); post=pre;
pre{1}=J; for j=1:M, pre{j+1}=mlpw.cascade(pre{j},seg{j}); end
post{M+1}=J; for j=M:-1:1, post{j}=mlpw.cascade(seg{j},post{j+1}); end
S=pre{end}; R.S_smm=S;
if strcmpi(method,'tmm')
 T=eye(2);
 for j=1:M
  a=seg{j}; v=a(2,1);
  if abs(v)<realmin, error('mlpw:TMMOverflow','TMM transmission underflow; use SMM.'); end
  T=T*[1/v -a(2,2)/v;a(1,1)/v a(1,2)-a(1,1)*a(2,2)/v];
 end
 assert(all(isfinite(T(:))),'TMM overflow; use SMM. No phase clipping is applied.');
 S=[T(2,1)/T(1,1),det(T)/T(1,1);1/T(1,1),-T(1,2)/T(1,1)]; R.transfer=T;
end
R.S=S; R.r=S(1,1); R.t=S(2,1); R.method=method;
R.C=C; R.q=q; R.k0=k0; R.z=[0 cumsum(d)];
R.Aleft=zeros(1,L); R.Bright=zeros(1,L);
% Unit left incidence; fields always reconstructed with stable S prefixes.
R.Aleft(1)=1; R.Bright(1)=R.S_smm(1,1);
R.Aleft(end)=R.S_smm(2,1);
for l=2:L-1
 P=pre{left(l)+1}; Q=post{left(l)+1};
 A=P(2,1)/(1-P(2,2)*Q(1,1));
 Qright=post{right(l)+1};
 R.Aleft(l)=A; R.Bright(l)=Qright(1,1)*A*exp(1i*q(l)*d(l-1));
end
% Physical port powers apply to lossless propagating entrance/exit only.
R.reflectance=NaN; R.transmittance=NaN; R.absorptance=NaN;
if isreal(e(1)) && isreal(e(end)) && real(e(1))>0 && real(e(end))>0 ...
        && isreal(C.kparallel) && isreal(q(1)) && isreal(q(end))
 R.reflectance=abs(R.r)^2; R.transmittance=real(q(end)/q(1))*abs(R.t)^2;
 R.absorptance=1-R.reflectance-R.transmittance;
end
end
