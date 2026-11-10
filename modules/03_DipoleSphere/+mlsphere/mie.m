function R=mie(C)
%MIE Standard electric a_l and magnetic b_l coefficients, passive sphere.
% C needs lambda/radius/epsilonSphere/nHost; position is unnecessary.
if ~isfield(C,'position'), C.position=[0;0;2*C.radius]; end
if ~isfield(C,'nmax')
 x=2*pi/C.lambda*C.nHost*C.radius; C.nmax=ceil(x+4*x^(1/3)+8);
end
S=mlsphere.model(C); B=mlsphere.modes(S);
R.l=B.l; R.a=B.a; R.b=B.b; R.nmax=S.nmax;
R.Qsca=2/S.x^2*sum((2*B.l+1).*(abs(B.a).^2+abs(B.b).^2));
R.Qext=2/S.x^2*sum((2*B.l+1).*real(B.a+B.b));
R.Qabs=R.Qext-R.Qsca;
end
