% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
% Claire Lin, 05/20/2018
clear all;close all;
load('data/cardiac_perf_R8.mat')
load('data/Xinf.mat')
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
opt.scaleL = 130/1.2775; % Otazo's stopping St(1) / b1 constant squared
opt.scaleS = 1/1.2775; % 1 / b1 constant squared
opt.muL=0.01;
opt.muS=0.01*opt.scaleS;
opt.Xinf = Xinf.perf;
%% AL-CG
d1 = 1/5; d2 = 1/5; %for AL-CG
[L_cg,S_cg,x_cg,cost_cg,time_cg,rankL_cg] = AL_CG(opt,'d1',d1,'d2',d2);
%% AL-2
d1 = 1/5; d2 = 1/50; %for AL-2
[L_al,S_al,xdiff_al,cost_al,time_al,rankL_al] = AL_2(opt,'d1',d1,'d2',d2);
%% prepare for PGM: param
param.E=getE(b1c,nt,'samp',kdata(:,:,:,1)~=0);
param.d=kdata;
param.T = getT(nx,ny,nt);
param.nite=10;
param.scaleL = 130/1.2775;
param.scaleS = 1/1.2775;
param.lambda_L=0.01;
param.lambda_S=0.01*param.scaleS;
param.Xinf = reshape(Xinf.perf,nx*ny,nt);
%% ISTA
[L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
%% FISTA
[L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1);
%% POGM
[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
%% Display: 4 frames
L = L_pogm;S = S_pogm;
LplusS=L+S;
LplusSd=LplusS(33:96,33:96,2);LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,8));LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,14));LplusSd=cat(2,LplusSd,LplusS(33:96,33:96,24));
Ld=L(33:96,33:96,2);Ld=cat(2,Ld,L(33:96,33:96,8));Ld=cat(2,Ld,L(33:96,33:96,14));Ld=cat(2,Ld,L(33:96,33:96,24));
Sd=S(33:96,33:96,2);Sd=cat(2,Sd,S(33:96,33:96,8));Sd=cat(2,Sd,S(33:96,33:96,14));Sd=cat(2,Sd,S(33:96,33:96,24));
figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,1]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,.03]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,1]);ylabel('S')