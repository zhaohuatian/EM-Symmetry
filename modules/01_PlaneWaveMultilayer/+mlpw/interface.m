function S=interface(e1,e2,q1,q2,pol)
%INTERFACE Unit electric-polarization amplitudes, exp(-i omega t).
if strcmpi(pol,'TE')
 r=(q1-q2)/(q1+q2); t=2*q1/(q1+q2); tr=2*q2/(q1+q2);
elseif strcmpi(pol,'TM')
 den=e2*q1+e1*q2; r=(e2*q1-e1*q2)/den;
 t=2*sqrt(e1)*sqrt(e2)*q1/den; tr=2*sqrt(e1)*sqrt(e2)*q2/den;
else, error('mlpw:Polarization','Use TE or TM.'); end
S=[r tr;t -r];
end
