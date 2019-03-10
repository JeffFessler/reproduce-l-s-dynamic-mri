function [L,S,x,cost,time,rankL] = PGM(param,varargin)
% [L,S,x,cost,time,rankL] = PGM(param,varargin)
% Claire Lin 05/20/2018
% use t=0.99 for ista, 0.5 for fista and pogm
%
% Inputs: 
% param.d: undersampled k-t data (nx,ny,nt,nc)
% param.E: data acquisition operator
% param.T: sparsifying transform
% param.lambda_L: nuclear-norm weight
% param.lambda_S: l1-norm weight
%
% Outputs: L, S, Xdifference, Cost, Time, rank of L

% arg
arg.fistaL = false; arg.fistaS = false;
arg.pogmL = false; arg.pogmS = false; 
arg.fista_restart = true;
arg.pogm_restart = true;
arg.L = param.E'*param.d; 
[nx,ny,nt]=size(arg.L);
arg.S = zeros(nx,ny,nt);  
arg.tscale = 1;
arg = vararg_pair(arg, varargin);
% initialize
Lpre = reshape(arg.L,[nx*ny,nt]); 
S=reshape(arg.S,[nx*ny,nt]); 
if arg.fistaL || arg.fistaS ||arg.pogmL || arg.pogmS
    t = 0.5*arg.tscale; 
else
    t = 0.99*arg.tscale; 
end
X = Lpre+S;
resLpre = param.E*reshape(X,[nx,ny,nt])-param.d;
M = X - t*reshape((param.E'*resLpre),[nx*ny,nt]);
YM = M; %for fista

if arg.fistaL 
    tL = t; YLpre = Lpre; ThetaLpre = 1; Thetapre = 1;
elseif arg.pogmL 
    tL = t; ULpre = Lpre; ZL = Lpre; ThetaLpre = 1; ZetaL = 1;
else
    tL = t;
end
if arg.fistaS
    tS = t; Spre = S; YS = S; ThetaSpre = 1;
elseif arg.pogmS
    tS = t; USpre = S; ZS = S; ThetaSpre = 1; ZetaS = 1;
else
    tS = t;
end
ite=0;
%print iter 0 cost
tmp2=param.T*reshape(S,[nx,ny,nt]);
x = zeros(1,param.nite+1); cost = x;
x(ite+1) = norm(col(X-param.Xinf),2);
cost(ite+1) = 0.5*norm(resLpre(:),2)^2+param.scaleL*param.lambda_L*sum(svd(reshape(Lpre,nx*ny,nt)))+param.lambda_S*norm(tmp2(:),1);
time = zeros(1,param.nite+1);
rankL = zeros(1,param.nite+1);rankL(1) = rank(Lpre);
fprintf(' ite: %d, xdiff: %f3, cost: %f3, rank of L: %d', ite,x(ite+1), cost(ite+1),rankL(ite+1)); 
	    
fprintf('\n ********** PGM: L+S reconstruction **********\n')
% iterations
for ite = 1:param.nite
	% L update
    tic;
    M0 = M;
    if arg.fistaL 
        [Ut,St,Vt]=svd(YM - YS,0);
        St=diag(SoftThresh(diag(St),param.scaleL*param.lambda_L*t));
        L=Ut*St*Vt';
    elseif arg.pogmL
        UL = M - S;
        if ite > param.nite - 1
            ThetaL = 0.5*(1+sqrt(1+8*ThetaLpre^2)); 
        else
            ThetaL = 0.5*(1+sqrt(1+4*ThetaLpre^2));
        end 
        c1 = (ThetaLpre-1)/ThetaL;
        c2 = ThetaLpre/ThetaL;
        c3 = -c1*tL/ZetaL;
        ZL = UL + c1*(UL-ULpre) + c2*(UL-Lpre) + c3*(Lpre-ZL);
        ZetaL = tL*(1+c1+c2); 
        [Ut,St,Vt]=svd(ZL,0);
        St=diag(SoftThresh(diag(St),param.scaleL*param.lambda_L*ZetaL)); 
        L=Ut*St*Vt';
        ULpre = UL;
        ThetaLpre = ThetaL;
    else
        [Ut,St,Vt]=svd(M - S,0);
        St=diag(SoftThresh(diag(St),param.scaleL*param.lambda_L*tL));
        L=Ut*St*Vt';
    end
    % S update
    if arg.fistaS
        S=reshape(param.T'*(SoftThresh(param.T*reshape(YM -  YLpre,[nx,ny,nt]),param.lambda_S*t)),[nx*ny,nt]);
    elseif arg.pogmS
        US = M - Lpre; 
        if ite > param.nite - 1
            ThetaS = 0.5*(1+sqrt(1+8*ThetaSpre^2)); 
        else
            ThetaS = 0.5*(1+sqrt(1+4*ThetaSpre^2));
        end
        ZS = US + (ThetaSpre-1)/ThetaS*(US-USpre) + ...
            ThetaSpre/ThetaS*(US-S) - ...
            (ThetaSpre-1)/ThetaS*tS/ZetaS*(S-ZS);
        ZetaS = tS*(1+(ThetaSpre-1)/ThetaS+ThetaSpre/ThetaS);
        S=reshape(param.T'*(SoftThresh(param.T*reshape(ZS,[nx,ny,nt]),param.lambda_S*ZetaS)),[nx*ny,nt]);
        USpre = US;
        ThetaSpre = ThetaS;
    else
        S=reshape(param.T'*(SoftThresh(param.T*reshape(M -  Lpre,[nx,ny,nt]),param.lambda_S*tS)),[nx*ny,nt]);
    end
    % data consistency
    X = L+S;
    resk=param.E*reshape(X,[nx,ny,nt])-param.d;
    M=X-t*reshape(param.E'*resk,[nx*ny,nt]); 
    % compute cost and update for printing and for restart
    tmp2=param.T*reshape(S,[nx,ny,nt]);
    x(ite+1) = norm(col(X-param.Xinf),2);
    cost(ite+1) = 0.5*norm(resk(:),2)^2+...
        param.scaleL*param.lambda_L*sum(svd(reshape(L,nx*ny,nt)))+param.lambda_S*norm(tmp2(:),1);
    % FISTA restart
    if arg.fistaL && arg.fistaS
        % update y
        if arg.fista_restart && cost(ite+1)>cost(ite) %function restart
            Theta = 1;
            YLpre = L;
            YS = S;
            printm('fista: restart at %d', ite)
        else
            Theta = (1+sqrt(1+4*Thetapre^2))/2;
            Beta = (Thetapre-1)/Theta;
            YLpre = L+Beta*(L-Lpre);
            YS = S+Beta*(S-Spre);
        end
        Spre = S;
        Thetapre = Theta;
        % data consistency
        reskY=param.E*reshape(YLpre+YS,[nx,ny,nt])-param.d;
        YM=YLpre+YS-t*reshape(param.E'*reskY,[nx*ny,nt]);
    end
	% POGM restart
    if arg.pogmL && arg.pogmS
        if arg.pogm_restart && cost(ite+1)>cost(ite) %function restart
            ThetaLpre = 1; ThetaSpre = 1;
            printm('pogm: restart at %d', ite)
        end
    end
    % update Lpre
	Lpre=L;
	% print cost function and solution update
    time(ite+1) = time(ite)+toc;
    rankL(ite+1) = rank(L);
    fprintf(' ite: %d, rank of L: %d, xdiff: %f3, cost: %f3\n',...
        ite,rankL(ite+1),x(ite+1),cost(ite+1));
end
L=reshape(L,nx,ny,nt);
S=reshape(S,nx,ny,nt);
end

% soft-thresholding function
function y=SoftThresh(x,p)
y=(abs(x)-p).*sign(x).*(abs(x)>p);
end    
