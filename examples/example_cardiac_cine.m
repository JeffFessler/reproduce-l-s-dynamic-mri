% L+S reconstruction of undersampled multicoil cardiac cine MRI
% Claire Lin, 05/20/2018
clear all;close all;
load('data/cardiac_cine_R6.mat')
load ('data/Xinf.mat')
%% normalize smap
[nx,ny,nt,nc]=size(kdata);
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
%% prepare for AL: opt
opt.d = kdata;
opt.smap = b1c;
opt.T=getT(nx,ny,nt);
opt.nite = 10;
opt.samp = opt.d(:,:,:,1)~=0;
[opt.F,opt.C] = getFS(opt.smap,nt,'samp',opt.samp);
opt.E=getE(b1c,nt,'samp',opt.samp);
% scalers to match Otazo's results
opt.scaleL = 175/1.2656; % Otazo's stopping St(1) / b1 constant squared
opt.scaleS = 1/1.2656; % 1 / b1 constant squared
opt.muL=0.01;
opt.muS=0.0025*opt.scaleS;
opt.Xinf = Xinf.cine;
%% AL-CG
d1 = 1/10; d2 = 1/20; %for AL-CG
[L_cg,S_alcg,x_alcg,cost_alcg,time_alcg,rankL_alcg] = AL_CG(opt,'d1',d1,'d2',d2);
%% AL-2
d1 = 1/10; d2 = 1/100;  %for AL-2
[L_al2,S_al2,xdiff_al2,cost_al2,time_al2,rankL_al2] = AL_2(opt,'d1',d1,'d2',d2);
%% prepare for PGM: param
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
param.d=kdata;
param.T = getT(nx,ny,nt);
param.nite=10;
param.scaleL = 175/1.2656;
param.scaleS = 1/1.2656;
param.lambda_L=0.01;
param.lambda_S=0.0025*param.scaleS;
param.Xinf = reshape(Xinf.cine,nx*ny,nt);
%% ISTA
[L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
%% FISTA
[L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1);
%% POGM
[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
%% Display: 4 frames
L = L_pogm;S = S_pogm;
LplusS=L+S;
% display 4 frames
LplusSd=LplusS(65:192,65:192,2);LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,8));LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,14));LplusSd=cat(2,LplusSd,LplusS(65:192,65:192,20));
Ld=L(65:192,65:192,2);Ld=cat(2,Ld,L(65:192,65:192,8));Ld=cat(2,Ld,L(65:192,65:192,14));Ld=cat(2,Ld,L(65:192,65:192,20));
Sd=S(65:192,65:192,2);Sd=cat(2,Sd,S(65:192,65:192,8));Sd=cat(2,Sd,S(65:192,65:192,14));Sd=cat(2,Sd,S(65:192,65:192,20));
figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,1]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,1]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,1]);ylabel('S')