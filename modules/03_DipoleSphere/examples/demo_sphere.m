function demo_sphere
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
C=struct('lambda',600e-9,'radius',80e-9,'epsilonSphere',4,'nHost',1, ...
 'position',[0;0;140e-9],'nmax',40); p=1e-29*[1;1i;.4];
gap=linspace(15,250,61)*1e-9; values=zeros(numel(gap),4);
for j=1:numel(gap)
 C.position=[0;0;C.radius+gap(j)]; C.nmax=max(40,ceil(log(1e-10)/(2*log(C.radius/norm(C.position)))));
 A=mlsphere.decay([0;0;1e-29],C); B=mlsphere.decay([1e-29;0;0],C);
 values(j,:)=[A.Ftotal A.Frad B.Ftotal B.Frad];
end
C.position=[0;0;140e-9]; C.nmax=40;
[t,f]=meshgrid(linspace(0,pi,101),linspace(0,2*pi,181)); F=mlsphere.farfield(p,C,t,f);
G=mlsphere.decay(p,C); fprintf('Sphere example Frad %.9g Ftotal %.9g\n',G.Frad,G.Ftotal);
fig=figure('Visible','on','Position',[100 100 1100 430]);
subplot(1,2,1); plot(gap*1e9,values(:,1),'r-',gap*1e9,values(:,3),'b-','LineWidth',1.4); grid on;
xlabel('Surface gap [nm]'); ylabel('Total decay / vacuum'); legend('Radial','Tangential'); title('Dielectric sphere: orientation dependence');
subplot(1,2,2); surf(F.pattern.*sin(t).*cos(f),F.pattern.*sin(t).*sin(f),F.pattern.*cos(t),F.dPdOmega,'EdgeColor','none');
axis equal; view(40,25); colorbar; xlabel('x'); ylabel('y'); zlabel('z'); title('Far-field shape; color dP/dOmega [W/sr]');
end
