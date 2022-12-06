% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
% Claire Lin, 05/20/2018
% 2022-12-05 Jeff Fessler

%% load data
if ~isvar('kdata')
    clear all; close all;
    addpath('../operators')
    addpath('../algorithms')

    load('../data/cardiac_perf_R8.mat')
    load('../data/Xinf.mat')
    Xinf_perf = Xinf.perf;
    clear Xinf
    [nx, ny, nt, nc] = size(kdata);
end

%% normalize smap
if ~isvar('b1c')
    ssos_raw = sqrt(sum(abs((b1)).^2,3)); % nearly uniform at 1.1303
    b1c = div0(b1, ssos_raw);
    im(b1c)
end

%% prepare for AL: opt
opt.d = kdata;
opt.smap = b1c;
opt.T = getT(nx,ny,nt);
opt.nite = 50;
opt.samp = kdata(:,:,:,1) ~= 0;
[opt.F, opt.C] = getFS(opt.smap, nt, 'samp', opt.samp);
opt.E = getE(b1c, nt, 'samp', opt.samp);

tmp = opt.E' * kdata; % adjoint (zero-filled) recon
%tmp = kdata;
%tmp(:,:,:,2:end) = 0; % 1st frame only
%tmp(:,:,2:end,:) = 0; % 1st coil only
%tmp = opt.E' * tmp; % adjoint (zero-filled) recon
im(tmp)
return

% scalars to match Otazo's results
opt.scaleL = 130/1.2775; % Otazo's stopping St(1) / b1 constant squared
opt.scaleS = 1/1.2775; % 1 / b1 constant squared
opt.muL = 0.01;
opt.muS = 0.01 * opt.scaleS;
opt.Xinf = Xinf_perf;

rerun = false;
%rerun = true

%% AL-CG
if ~isvar('L_cg') % || rerun
    d1 = 1/5; d2 = 1/5; % for AL-CG
    [L_cg,S_cg,x_cg,cost_cg,time_cg,rankL_cg] = AL_CG(opt,'d1',d1,'d2',d2);
end

%% AL-2
if ~isvar('L_al') % || rerun
    d1 = 1/5; d2 = 1/50; % for AL-2
    [L_al,S_al,xdiff_al,cost_al,time_al,rankL_al] = AL_2(opt,'d1',d1,'d2',d2);
end


%% prepare for PGM: param
%param.E = getE(b1c, nt, 'samp', opt.samp); % todo cut
param.E = opt.E;
param.d = kdata;
param.T = getT(nx,ny,nt);
param.nite = opt.nite;
param.scaleL = 130/1.2775;
param.scaleS = 1/1.2775;
param.lambda_L = 0.01;
param.lambda_S = 0.01*param.scaleS;
param.Xinf = reshape(Xinf_perf, nx*ny, nt);
%return

%% ISTA
if ~isvar('L_ista') || rerun
    [L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
end

%% FISTA
if ~isvar('L_fista') || rerun
    [L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = ...
    PGM(param,'fistaL',1,'fistaS',1);
end

%% POGM
if ~isvar('L_pogm') || rerun
    [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = ...
    PGM(param,'pogmS',1,'pogmL',1);
end

%% Display: 4 frames
ix = 33:96; iy = 33:96; % roi
it = [2 8 14 24]; % frames
mycat = @(x) cat(2, x(:,:,1), x(:,:,2), x(:,:,3), x(:,:,4));
L = L_pogm; S = S_pogm;
LplusS = L + S;
LplusSd = mycat(LplusS(ix,iy,it));
Ld = mycat(L(ix,iy,it));
Sd = mycat(S(ix,iy,it));
figure(3)
subplot(3,1,1), imshow(abs(LplusSd), [0,1]); ylabel('L+S')
subplot(3,1,2), imshow(abs(Ld), [0,.03]); ylabel('L')
subplot(3,1,3), imshow(abs(Sd), [0,1]); ylabel('S')

%% Plot cost function (cf. Fig. 1)
% (In the paper it is relative to the cost at infinity,
% but we cannot compute that from Xinf alone; need Linf and Sinf.)
niter = param.nite;
iter = 0:niter;

fun = @(x) x / cost_cg(1);
figure(2)
semilogy(...
 iter, fun(cost_cg), 'g-o', ...
 iter, fun(cost_ista), 'b-x', ...
 iter, fun(cost_fista), 'k-*', ...
 iter, fun(cost_al), 'm-^', ...
 iter, fun(cost_pogm), 'r-s' ...
)
xlabel('iteration'); ylabel('cost / cost_0')

figure(1)
semilogy(...
 time_cg,    fun(cost_cg),    'g-o', ...
 time_ista,  fun(cost_ista),  'b-x', ...
 time_fista, fun(cost_fista), 'k-*', ...
 time_al,    fun(cost_al),    'm-^', ...
 time_pogm,  fun(cost_pogm),  'r-s' ...
)
legend('AL-CG', 'ISTA', 'FISTA', 'AL-2', 'POGM'), grid
xlabel('wall time (s)'); ylabel('cost / cost_0')
xlim([0,30])
