% L+S reconstruction of undersampled multicoil PINCAT phantom
% Claire Lin, 06/05/2018
clear all;close all;
load('data/aperiodic_pincat.mat')
load('data/Xinf.mat')
tmp = max(new(:));
Xtrue = div0(new,tmp);
[nx,ny,nt] = size(Xtrue);
im(Xtrue)
%% simulate coil maps
% 4 rings of 8 coils
nc = 8;
nring = 4;
b1 = ir_mri_sensemap_sim('chat', 0, 'nx', nx, 'ny', 1, 'nz', ny,'dx', 1.5, 'dz', 1.5, ...
	'rcoil', 120, 'nring', nring, 'ncoil', nc*nring, 'coil_distance', 1.2);
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
im(b1c)
%% coil compression to 8 coils
[smap0, ~] = ir_mri_coil_compress(b1c, 'ncoil', nc);
tmp = sqrt(sum(abs((smap0)).^2,3));
smap = div0(smap0,tmp);
im(smap)
%% sampling mask
line = 24;
samp = strucrand(nx,ny,nt,line);
im(samp)
%% data
E=getE(smap,nt,'samp',samp);
d = E*Xtrue;
% add noise
rng(0)
dn = randn(size(d)) + 1i * randn(size(d));
param.snr_db = 46;
param.scale_noise = norm(d(:)) / norm(dn(:)) / 10.^(param.snr_db / 20);
param.d = d + param.scale_noise * dn;
printm('data rmse = %g, snr = %g', rms(param.d(:)-d(:)), ...
    20*log10(norm(d(:)) / norm(param.d(:)-d(:))))
opt.d = param.d;
% prepare for regularization scaling
L = E'*param.d; 
res = E*L-param.d;
[~,St,~]=svd(reshape(L,[nx*ny,nt])-reshape(E'*res,[nx*ny,nt]),0);
semilogy(diag(St),'ok')
%% prepare for AL: opt
opt.smap = smap;
opt.T=getT(nx,ny,nt);
opt.nite = 10;
opt.samp = samp;
[opt.F,opt.C] = getFS(opt.smap,nt,'samp',opt.samp);
opt.E=getE(smap,nt,'samp',opt.samp);
opt.scaleL = St(1);
opt.scaleS = 1/1.887; %1 / b1 constant squared in middle of image
opt.muL=0.01;
opt.muS=0.05*opt.scaleS;
opt.Xinf = Xinf.pincat;
%% AL-CG
d1 = 1/2; d2 = 1/2; %for AL-CG
[L_cg,S_cg,x_cg,cost_cg,time_cg,rankL_cg] = AL_CG(opt,'d1',d1,'d2',d2); 
%% AL-2
d1 = 1/3; d2 = 1/10; %for AL-2
[L_al,S_al,xdiff_al,cost_al,time_al,rankL_al] = AL_2(opt,'d1',d1,'d2',d2);
%% prepare for PGM: param
param.E = getE(smap,nt,'samp',samp);
param.T = getT(nx,ny,nt);
param.nite = 10;
param.scaleL = St(1);
param.scaleS = 1/1.887; %1 / b1 constant squared in middle of image
param.lambda_L=0.01;
param.lambda_S=0.05*param.scaleS;
param.Xinf = reshape(Xinf.pincat,nx*ny,nt);
%% ISTA
[L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
%% FISTA
[L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1);
%% POGM
[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
%% Display: 1 frame
figure; tpick = 21; L = L_pogm;S = S_pogm;
LplusS=L+S;
subplot(1,5,1)
imshow(abs(Xtrue(:,:,tpick)),[0,1])
title({'$X_{true}$'},'Interpreter','latex')
subplot(1,5,2)
imshow(abs(LplusS(:,:,tpick)),[0,1])
title({'$\tilde{X}$'},'Interpreter','latex')
subplot(1,5,3)
imshow(abs(L(:,:,tpick)),[0,1])
title({'$L_\infty$'},'Interpreter','latex')
subplot(1,5,4)
imshow(abs(S(:,:,tpick)),[0,1])
title({'$S_\infty$'},'Interpreter','latex')
subplot(1,5,5)
imagesc(abs(LplusS(:,:,21))-abs(Xtrue(:,:,tpick)),[0,0.2])
title({'$|X_{true}-\tilde{X}|$'},'Interpreter','latex')
axis off;axis square;