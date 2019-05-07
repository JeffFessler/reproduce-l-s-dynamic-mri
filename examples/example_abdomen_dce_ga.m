% L+S reconstruction of undersampled multicoil abdominal dce MRI
% Claire Lin, 05/07/2019
%clear all;close all;
load('data/abdomen_dce_ga.mat')
load ('data/Xinf.mat')
%% number of spokes to be used per frame (Fibonacci number)
nspokes=21;
[nx,ny,nc]=size(b1);
[nr,ntviews,nc]=size(kdata);
% number of frames
nt=floor(ntviews/nspokes);
% crop the data according to the number of spokes per frame
kdata=kdata(:,1:nt*nspokes,:);
k=k(:,1:nt*nspokes);
w=w(:,1:nt*nspokes);
% sort the data into a time-series of undersampled images
for ii=1:nt
    kdatau(:,:,:,ii)=kdata(:,(ii-1)*nspokes+1:ii*nspokes,:);
    ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end
%% normalize smap
tmp = sqrt(sum(abs((b1)).^2,3));
b1c = div0(b1,tmp);
%% prepare for system operator
ksp = reshape(ku,[],nt); 
ksp3(:,1,:) = real(ksp);
ksp3(:,2,:) = imag(ksp);
om3 = ksp3*2*pi;
wi3 = reshape(wu,[],nt);
M = size(ksp,1);
%% prepare for PGM: param
param.E = getEnufft(b1c,nt,'M',M,'ksp',ksp3,'om',om3,'wi',wi3,'donufft',1);
param.d=permute(reshape(kdatau,[],nc,nt),[1 3 2]);
recon_nufft=param.E'*param.d;
tscale = 1.3002; %power iteration on E'E
param.scaleL = 0.0274/1.099;
param.scaleS = 1/1.099;
param.lambda_L=0.025;
param.lambda_S=0.5*max(abs(recon_nufft(:)))*param.scaleS;
param.nite=10;
param.tol=0.0025;
param.T= getT(nx,ny,nt);
param.Xinf = reshape(Xinf.abd,nx*ny,nt);
%% ISTA
[L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param,'tscale',tscale);
%% FISTA
[L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1,'tscale',tscale);
%% POGM
[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1,'tscale',tscale);
%% Display: 4 frames
L = L_pogm;S = S_pogm;
L=flipdim(L,1);
S=flipdim(S,1);
LplusS=L+S;
% display 4 frames
LplusSd=LplusS(65:336,:,1);LplusSd=cat(2,LplusSd,LplusS(65:336,:,9));LplusSd=cat(2,LplusSd,LplusS(65:336,:,16));LplusSd=cat(2,LplusSd,LplusS(65:336,:,25));
Ld=L(65:336,:,1);Ld=cat(2,Ld,L(65:336,:,9));Ld=cat(2,Ld,L(65:336,:,16));Ld=cat(2,Ld,L(65:336,:,25));
Sd=S(65:336,:,1);Sd=cat(2,Sd,S(65:336,:,9));Sd=cat(2,Sd,S(65:336,:,16));Sd=cat(2,Sd,S(65:336,:,25));
figure;
subplot(3,1,1),imshow(abs(LplusSd),[0,5e-4]);ylabel('L+S')
subplot(3,1,2),imshow(abs(Ld),[0,5e-4]);ylabel('L')
subplot(3,1,3),imshow(abs(Sd),[0,5e-4]);ylabel('S')
