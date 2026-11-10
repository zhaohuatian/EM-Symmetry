function demo_planewave
%DEMO_PLANEWAVE Thin-film spectrum and a physical Cartesian standing wave.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
C=struct('lambda',600e-9,'epsilon',[1 2.1^2 1.45^2 1], ...
 'thickness',[90 160]*1e-9,'kparallel',0,'polarization','TE');
wl=linspace(400,900,201)*1e-9; values=zeros(numel(wl),3);
for j=1:numel(wl)
 C.lambda=wl(j); R=mlpw.solve(C); values(j,:)=[R.reflectance R.transmittance R.absorptance];
end
C.lambda=600e-9; R=mlpw.solve(C); F=mlpw.field(R,linspace(-200e-9,sum(C.thickness)+200e-9,801));
fig=figure('Visible','on','Position',[100 100 1100 430]);
subplot(1,2,1); plot(wl*1e9,values,'LineWidth',1.3); grid on; xlabel('Vacuum wavelength [nm]'); ylabel('Power fraction'); legend('R','T','A'); title('Lossless multilayer spectrum');
subplot(1,2,2); plot(F.z*1e9,sum(abs(F.E).^2,2),'LineWidth',1.3); grid on; xlabel('z [nm]'); ylabel('|E|^2 / |E_{inc}|^2'); title('Physical electric field at 600 nm');
for zz=R.z, xline(zz*1e9,'--'); end
end
