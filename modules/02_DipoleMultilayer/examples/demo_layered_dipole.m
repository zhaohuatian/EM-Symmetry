function demo_layered_dipole
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
D=struct('WL0',655e-9,'Eplist',[1.4564^2 1 1],'dl',[0 100e-9], ...
 'nDn',1,'POSD',17e-9); p=1e-29*[1;1i;.4];
R=mldipole.radiative(p,D,120); G=mldipole.decay(p,D);
fprintf('Layered example Frad %.9g Ftotal %.9g\n',R.Frad,G.Ftotal);
[t,f]=meshgrid(linspace(0,pi/2-1e-5,101),linspace(0,2*pi,181));
fig=figure('Visible','on','Position',[100 100 1100 430]);
for side=1:2
 sgn=3-2*side; F=mldipole.farfield(p,D,t,f,sgn);
 subplot(1,2,side); surf(sin(t).*cos(f),sin(t).*sin(f),F.dPdOmega,'EdgeColor','none');
 view(2); axis equal tight; colorbar; xlabel('sin(theta) cos(phi)'); ylabel('sin(theta) sin(phi)');
 if sgn==1, title('Upper dP/dOmega [W/sr]'); else, title('Lower dP/dOmega [W/sr]'); end
end
end
