function [F,S] = getFS(sense_maps,Nt,varargin)
% Claire Lin 05/20/2018
% outputs: F = Omega, S = QC 
arg.Nx = size(sense_maps,1);
arg.Ny = size(sense_maps,2);
arg.Nt = Nt;
arg.Nc = size(sense_maps,3);
arg.samp = []; 
arg.smaps = sense_maps;
arg.fullmask = true(arg.Nx, arg.Ny, arg.Nt, arg.Nc);
arg.submask = arg.fullmask(:,:,:,1); %assume mask same across coils
arg.N = [arg.Nx,arg.Ny];
arg.fov = 22;
arg.basis = {'dirac'};
arg.ksp = mri_trajectory('cartesian', {}, arg.N, arg.fov);
arg.dim = [size(arg.ksp,1) sum(arg.fullmask(:))];
arg = vararg_pair(arg, varargin);

if (arg.Nc > 1)
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt arg.Nc], ...
                    'arg', arg, 'odim', [arg.Nx arg.Ny arg.Nt arg.Nc],... 
                    'forw', @F_forw, 'back', @F_back);
else 
                F = fatrix2('idim', [arg.Nx arg.Ny arg.Nt], ...
                    'arg', arg,'odim', [arg.Nx arg.Ny arg.Nt],...
                    'forw', @F_forw, 'back', @F_back);
end
S = fatrix2('idim', [arg.Nx arg.Ny arg.Nt] ,'arg',arg,...
    'odim', [arg.Nx arg.Ny arg.Nt arg.Nc],...
    'forw', @S_forw,'back', @S_back);
end

% y = A * x
function S = F_forw(arg,s) 
if ~isempty(arg.samp)
    for ch=1:arg.Nc
        S(:,:,:,ch)=s(:,:,:,ch).*arg.samp;
    end
end
end

function s = S_forw(arg,x)
for tt=1:arg.Nt
    for ch=1:arg.Nc
        s(:,:,tt,ch)=x(:,:,tt).*arg.smaps(:,:,ch); 
    end
end
s=fft2c_mri(s);
end

% x = A' * y
function s = F_back(arg,S)
if ~isempty(arg.samp)
    for ch=1:arg.Nc
        s(:,:,:,ch)=S(:,:,:,ch).*arg.samp;
    end
end
end

function x = S_back(arg,s)
s = ifft2c_mri(s);
for tt=1:arg.Nt
    x(:,:,tt)=sum(squeeze(s(:,:,tt,:)).*conj(arg.smaps),3);
end
end