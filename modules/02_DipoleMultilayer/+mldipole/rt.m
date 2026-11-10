function [Ru,Rd,Tu,Td,q]=rt(S,rho)
%RT Adapt the dipole geometry to module 1; no duplicate Fresnel recursion.
% mlpw.solve returns electric s/p amplitudes at stack interface planes.
q=sqrt(S.k.^2-rho.^2); q(imag(q)<0)=-q(imag(q)<0);
[Ru,Tu]=side(1); [Rd,Td]=side(-1);
    function [R,T]=side(direction)
        exterior=1+(direction>0)*(numel(S.ep)-1);
        ids=S.j:direction:exterior;
        R=[0 0]; T=[1 1];
        if numel(ids)==1, return; end
        finite=ids(2:end-1);
        C=struct('lambda',S.WL0,'epsilon',S.ep(ids), ...
          'thickness',S.z(finite)-S.z(finite-1), ...
          'kparallel',rho,'polarization','TE');
        for polarization=1:2
            if polarization==2, C.polarization='TM'; end
            response=mlpw.solve(C,'smm');
            R(polarization)=response.r; T(polarization)=response.t;
        end
        if direction>0, distance=S.z(S.j)-S.POSD;
        else, distance=S.POSD-S.z(S.j-1); end
        phase=exp(1i*q(S.j)*distance);
        R=R*phase^2; T=T*phase;
    end
end
