function R=decay(p,D,opt)
% Total electromagnetic decay from regularized local Green tensor.
% Contour below real k_parallel captures lossless guided poles causally.
if nargin<3, opt=struct; end
S=mldipole.model(D); p=p(:);
assert(numel(p)==3 && all(isfinite(p)) && norm(p)>0,'Invalid p.');
if ~isfield(opt,'eta'), opt.eta=0.02; end
if ~isfield(opt,'RelTol'), opt.RelTol=1e-7; end
if ~isfield(opt,'AbsTol'), opt.AbsTol=1e-9; end
dist=abs(S.z-S.POSD);
if ~isfield(opt,'umax'), opt.umax=max([10,2*max(abs(sqrt(S.ep))),30/(S.k0*min([dist inf]))]); end
assert(opt.eta>0 && opt.umax>max(abs(sqrt(S.ep))),'Invalid contour.');
gx=integral(@(t) fun(t,1),0,opt.umax,'RelTol',opt.RelTol,'AbsTol',opt.AbsTol);
gz=integral(@(t) fun(t,2),0,opt.umax,'RelTol',opt.RelTol,'AbsTol',opt.AbsTol);
% Complete finite contour with vertical endpoint segment.
uend=opt.umax-1i*opt.eta*opt.umax/(1+opt.umax);
gx=gx+integral(@(v) endpoint(v,1),0,1,'AbsTol',opt.AbsTol,'RelTol',opt.RelTol);
gz=gz+integral(@(v) endpoint(v,2),0,1,'AbsTol',opt.AbsTol,'RelTol',opt.RelTol);
R.Gsc=S.k0*diag([gx gx gz]);
R.ImG=S.nh*S.k0/(6*pi)*eye(3)+imag(R.Gsc);
R.Ftotal=6*pi/S.k0*real(p'*R.ImG*p)/sum(abs(p).^2);
R.P0=S.omega^4*sum(abs(p).^2)/(12*pi*S.eps0*S.c^3);
R.Ptotal=R.P0*R.Ftotal;
R.Esc=S.mu0*S.omega^2*R.Gsc*p;
R.Gamma0=S.omega^3*sum(abs(p).^2)/(3*pi*S.eps0*1.054571817e-34*S.c^3);
R.Gamma=R.Gamma0*R.Ftotal; R.options=opt;
    function y=fun(t,c)
        u=t-1i*opt.eta*t./(1+t); du=1-1i*opt.eta./(1+t).^2;
        y=arrayfun(@(a) kernel(a,c),u).*du;
    end
    function y=endpoint(v,c)
        du=opt.umax-uend;
        y=arrayfun(@(a) kernel(a,c),uend+v*du)*du;
    end
    function g=kernel(u,c)
        [a,b,~,~,q]=mldipole.rt(S,S.k0*u); w=q(S.j)/S.k0;
        den=1-a.*b;
        plus=(a+b+2*a.*b)./den;
        if c==1
            minus=(a+b-2*a.*b)./den;
            g=1i/(8*pi)*u/w*(plus(1)-w^2/S.ep(S.j)*minus(2));
        else
            g=1i/(4*pi)*u^3/(w*S.ep(S.j))*plus(2);
        end
    end
end
