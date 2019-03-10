function [L,S,x,cost,time,rankL] = AL_2(opt,varargin)
% [L,S,x,cost,time,rankL] = AL_2(opt,varargin)
% Claire Lin 05/20/2018
% 4 variables, efficient implementation
% split: F - sampling mask & C - Fourier and coil sensitivity
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
samp = opt.samp;
F = opt.F;
C = opt.C;
arg.L = C'*(F'*opt.d); 
arg.S = zeros(nx,ny,nt);
arg.d1 = 1/10; arg.d2 = 1/100; 
arg = vararg_pair(arg, varargin);
% initialization
L = arg.L; S = arg.S; X = L+S; 
V1 = zeros(size(opt.d)); 
Flhs = 1./(samp+arg.d1);
ceig = sum(abs(opt.smap).^2,3);
V2 = zeros(size(L));
res = F*(C*X)-opt.d;
% print iter 0 cost
ite=0;
tmp2=opt.T*S;
x = zeros(1,opt.nite+1); cost = x;
x(ite+1) = norm(col(X-opt.Xinf),2);
cost(ite+1) = 0.5*norm(res(:),2)^2+opt.scaleL*opt.muL*sum(svd(reshape(L,nx*ny,nt)))+opt.muS*norm(tmp2(:),1);
time = zeros(1,opt.nite+1);
rankL = zeros(1,opt.nite+1);rankL(1) = rank(reshape(L,nx*ny,nt));
fprintf(' ite: %d, xdiff: %f3, cost: %f3, rank of L: %d', ite,x(ite+1),cost(ite+1),rankL(1)); 
fprintf('\n ********** AL-2: L+S reconstruction **********\n')
% loop through updates of L and S
for ite = 1:opt.nite
    tic;
    % updates
    % 1:Z 
    rhs = F'*opt.d+arg.d1*(C*X-V1);
    Z = bsxfun(@times,Flhs,rhs); 
    % 2:X
    rhs = sum(C'*(Z+V1),4) + arg.d2/arg.d1*(L+S-V2); 
    X = bsxfun(@rdivide,rhs,ceig+arg.d2/arg.d1);
    % 3:L
    [Ut,St,Vt] = svd(reshape(X-S+V2,[nx*ny,nt]),0);
    St = diag(SoftThresh(diag(St),opt.scaleL*opt.muL/(arg.d2)));
    L = reshape(Ut*St*Vt',[nx,ny,nt]);
    % 4:S
    S = opt.T'*SoftThresh(opt.T*(X-L+V2),opt.muS/(arg.d2));
    % 5:V1,V2
    V1 = V1 + (Z-C*X);
    V2 = V2 + (X-L-S);
    % print cost 
    time(ite+1) = time(ite)+toc;
    LpS=L+S;
    res = F*(C*(LpS))-opt.d;
    tmp2=opt.T*S;
    rankL(ite+1) = rank(reshape(L,[nx*ny,nt]));
    x(ite+1) = norm(col(LpS-opt.Xinf),2);
    cost(ite+1) = 0.5*norm(res(:),2)^2+opt.scaleL*opt.muL*sum(svd(reshape(L,nx*ny,nt)))+opt.muS*norm(tmp2(:),1);
    fprintf(' ite: %d , rank of L: %d, xdiff: %f3, cost: %f3\n', ite,rankL(ite+1),x(ite+1),cost(ite+1));
end
end

function y = SoftThresh(x,p)
    y = (abs(x)-p).*sign(x).*(abs(x)>p);
end   
