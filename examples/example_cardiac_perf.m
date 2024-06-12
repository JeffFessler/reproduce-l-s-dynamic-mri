% L+S reconstruction of undersampled multicoil cardiac perfusion MRI
% Claire Lin, 05/20/2018
% 2022-12-05 Jeff Fessler

%% load data (download data from web if needed)
if ~isvar('kdata')
    clear all; close all;
    addpath('../operators')
    addpath('../algorithms')

    url = 'https://web.eecs.umich.edu/~fessler/irt/reproduce/19/lin-19-edp/data/';
    filename = 'cardiac_perf_R8.mat';
    filepath = ['../data/' filename];
    if ~exist(filepath, 'file')
          fileurl = [url filename];
          websave(filepath, fileurl);
    end
    load(filepath)

    filename = 'Xinf.mat';
    filepath = ['../data/' filename];
    if ~exist(filepath, 'file')
          fileurl = [url filename];
          websave(filepath, fileurl);
    end
    load(filepath)

    Xinf = Xinf.perf;
    clear filepath filename fileurl url
    [nx, ny, nt, nc] = size(kdata);
end


%% normalize smap
if ~isvar('b1c')
    ssos_raw = sqrt(sum(abs(b1).^2, 3)); % nearly uniform at 1.1303
    b1c = div0(b1, ssos_raw);
    clear b1 % make sure to use the normalized one after this
    jim(b1c)
%   extrema(@abs, sqrt(sum(abs(b1c).^2, 3))) % verify
end


%% prepare for AL: opt
opt.d = kdata;
opt.smap = b1c;
opt.T = getT(nx,ny,nt);
opt.nite = 10;
opt.samp = kdata(:,:,:,1) ~= 0;
[opt.F, opt.C] = getFS(opt.smap, nt, 'samp', opt.samp);
opt.E = getE(b1c, nt, 'samp', opt.samp);

%tmp = fftshift(opt.T * Xinf, 3); % for julia check
%tmp = opt.T' * (opt.T * Xinf); extrema(@abs, tmp - Xinf)
%tmp = opt.E * Xinf;
%tmp = permute(tmp, [1 2 4 3]); % time last for julia check
%tmp = single(opt.E' * kdata); % adjoint (zero-filled) recon
%tmp = opt.E * Xinf;
%tmp = tmp(:)'*kdata(:) / sum(abs2(tmp(:)))
%jim(tmp)
%save tmp.mat tmp

% scalars to match Otazo's results
opt.scaleL = 130/1.2775; % Otazo's stopping St(1) / b1 constant squared
opt.scaleS = 1/1.2775; % 1 / b1 constant squared
opt.muL = 0.01; % AL parameters
opt.muS = 0.01 * opt.scaleS;
opt.Xinf = Xinf;

rerun = false;
%rerun = true

%% AL-CG
if ~isvar('L_cg') % || rerun
    d1 = 1/5; d2 = 1/5; % for AL-CG
    [L_cg, S_cg, x_cg, cost_cg, time_cg, rankL_cg] ...
        = AL_CG(opt, 'd1', d1, 'd2', d2);
end

%% AL-2
if ~isvar('L_al') % || rerun
    d1 = 1/5; d2 = 1/50; % for AL-2
    [L_al, S_al, xdiff_al, cost_al, time_al, rankL_al] = ...
        AL_2(opt, 'd1', d1, 'd2', d2);
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
param.Xinf = reshape(Xinf, nx*ny, nt);

%% ISTA
if ~isvar('L_ista') || rerun
    [L_ista, S_ista, xdiff_ista, cost_ista, time_ista, rankL_ista] = PGM(param);
end

%% FISTA
if ~isvar('L_fista') || rerun
    [L_fista, S_fista, xdiff_fista, cost_fista, time_fista, rankL_fista] = ...
    PGM(param, 'fistaL', 1, 'fistaS', 1);
end

%% POGM
if ~isvar('L_pogm') || rerun
    [L_pogm, S_pogm, xdiff_pogm, cost_pogm, time_pogm, rankL_pogm] = ...
    PGM(param, 'pogmS', 1, 'pogmL', 1);
end

%% Display: 4 frames
ix = 33:96; iy = 33:96; % roi
it = [2 8 14 24]; % frames
mycat = @(x) cat(2, x(:,:,1), x(:,:,2), x(:,:,3), x(:,:,4));
myroi = @(x) mycat(x(ix,iy,it));
L = L_pogm; S = S_pogm;
LplusS = L + S;
LplusSd = myroi(LplusS);
Ld = myroi(L);
Sd = myroi(S);
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
