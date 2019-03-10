function [L,S,x,cost,time,rankL] = AL_CG(opt,varargin)
% [L,S,x,cost,time,rankL] = AL_CG(opt,varargin)
% Claire Lin 05/20/2018
% 4 variables, AL with CG
% Input: struct opt
% opt.d: undersampled k-t data (nx,ny,nt,nc)
% opt.smap: smap (nx,ny,nc)
% opt.T: sparsifying transform (temporal: 3rd dim)
% opt.muL: nuclear-norm reg param
% opt.muS: l1-norm reg param
% opt.nite: number of iter
% Outputs: L, S, Xdifference, Cost, Time, rank of L
%
%parameters
[nx,ny,nt,nc]=size(opt.d);
Etd = opt.E'*opt.d;
arg.L = Etd; 
arg.S = zeros(nx,ny,nt);
arg.d1 = 1/10; arg.d2 = 1/10; 
arg.iterL = 3; arg.iterS = 3;
arg = vararg_pair(arg, varargin);
% precompute
EtE = opt.E'*opt.E;
d1eye = Gmatrix(diag_sp(arg.d1*ones(nx*ny*nt,1)),'idim',[nx,ny,nt],'odim',[nx,ny,nt]);
d2eye = Gmatrix(diag_sp(arg.d2*ones(nx*ny*nt,1)),'idim',[nx,ny,nt],'odim',[nx,ny,nt]);
E1 = EtE+d1eye;
E2 = EtE+d2eye;
% initialization
L = arg.L; S = arg.S; 
V1 = zeros(size(L)); V2 = V1;
% reg param
res = opt.E*(L+S)-opt.d;
% print iter 0 cost
ite=0;
tmp2=opt.T*S;
x = zeros(1,opt.nite+1); cost = x;
x(ite+1) = norm(col(L+S-opt.Xinf),2);
cost(ite+1) = 0.5*norm(res(:),2)^2+opt.scaleL*opt.muL*sum(svd(reshape(L,nx*ny,nt)))+opt.muS*norm(tmp2(:),1);
rankL = zeros(1,opt.nite+1);rankL(1) = rank(reshape(L,nx*ny,nt));
fprintf(' ite: %d, xdiff: %f3, cost: %f3, rank of L: %d', ite,x(ite+1), cost(ite+1),rankL(1)); 
time = zeros(1,opt.nite+1);
fprintf('\n ********** AL-CG: L+S reconstruction **********\n')
for ite = 1:opt.nite
    tic;
    % updates
    % 1:P
    [Ut,St,Vt] = svd(reshape(L+V1,[nx*ny,nt]),0);
    St = diag(SoftThresh(diag(St),opt.scaleL*opt.muL/(arg.d2)));
    P = reshape(Ut*St*Vt',[nx,ny,nt]);
    % 2:Q
    Q = SoftThresh(opt.T*S+V2,opt.muS/(arg.d2));
    % 3:L
    L = cgs(@(z) E1*z,col(Etd-opt.E'*(opt.E*S)+arg.d1*(P-V1)),1e-5,arg.iterL,[],[],L(:));
    L = reshape(L,[nx,ny,nt]);
    % 4:S
    S = cgs(@(z) E2*z,col(Etd-opt.E'*(opt.E*L)+arg.d2*opt.T'*(Q-V2)),1e-5,arg.iterS,[],[],S(:));
    S = reshape(S,[nx,ny,nt]);
    % 5:V1,V2
    V1 = V1 + (L-P);
    V2 = V2 + (opt.T*S-Q);
    % print cost 
    time(ite+1) = time(ite)+toc;
    LpS = L+S;
    res = opt.E*LpS-opt.d;
    tmp2=opt.T*S;
    rankL(ite+1) = rank(reshape(L,[nx*ny,nt]));
    x(ite+1) = norm(col(LpS-opt.Xinf),2);
    cost(ite+1) = 0.5*norm(res(:),2)^2+opt.scaleL*opt.muL*sum(svd(reshape(L,nx*ny,nt)))+opt.muS*norm(tmp2(:),1);
    fprintf(' ite: %d , rank of L: %d, xdiff: %f3, cost: %f3\n', ite,rankL(ite+1),x(ite+1),cost(ite+1));
end

function y=SoftThresh(x,p)
y=(abs(x)-p).*sign(x).*(abs(x)>p);
end   

end