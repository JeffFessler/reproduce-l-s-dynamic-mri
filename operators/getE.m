function E = getE(sense_maps,Nt,varargin)
% Claire Lin 05/20/2018
arg.Nx = size(sense_maps,1);
arg.Ny = size(sense_maps,2);
arg.Nt = Nt;
arg.Nc = size(sense_maps,3);
arg.mask = true(arg.Nx, arg.Ny, arg.Nt);
arg.smaps = sense_maps;
arg.samp = []; 
arg.N = [arg.Nx,arg.Ny];
arg.fov = [22,22];
arg.ksp = mri_trajectory('cartesian', {}, arg.N, arg.fov);
arg.dim = [size(arg.ksp,1) sum(arg.mask(:))];
arg = vararg_pair(arg, varargin);
   

if (arg.Nc > 1)
    E = fatrix2('arg',arg,'imask', arg.mask, ...
                         'odim', [arg.Nx arg.Ny arg.Nt arg.Nc], 'forw', ...
                        @E_forw, 'back', @E_back);
else 
    E = fatrix2('arg',arg ,'imask', arg.mask,...
                         'odim', [arg.Nx arg.Ny arg.Nt], 'forw', ...
                        @E_forw, 'back', @E_back);
end
end

% y = A * x
function S = E_forw(arg,x)
s = zeros(arg.Nx,arg.Ny,arg.Nt,arg.Nc); 
% for tt=1:arg.Nt
%     for ch=1:arg.Nc
%         s(:,:,tt,ch)=x(:,:,tt).*arg.smaps(:,:,ch); 
%     end
% end
s = bsxfun(@times,x,reshape(arg.smaps,[arg.Nx,arg.Ny,1,arg.Nc]));
S=fft2c_mri(s);
% if ~isempty(arg.samp)
%     for ch=1:arg.Nc
%         S(:,:,:,ch)=S(:,:,:,ch).*arg.samp;
%     end
% end
if ~isempty(arg.samp) %cl: samp mask only when cartesian samp
    S = bsxfun(@times,S,arg.samp);
end
end

% x = A' * y
function x = E_back(arg,S)
x = zeros(arg.Nx,arg.Ny,arg.Nt);
if ~isempty(arg.samp)
%     for ch=1:arg.Nc
%         S(:,:,:,ch)=S(:,:,:,ch).*arg.samp;
%     end
    S = bsxfun(@times,S,arg.samp);
end
s = ifft2c_mri(S);
% for tt=1:arg.Nt
%     x(:,:,tt)=sum(squeeze(s(:,:,tt,:)).*conj(arg.smaps),3);
% end
x = sum(bsxfun(@times,s,reshape(conj(arg.smaps),[arg.Nx,arg.Ny,1,arg.Nc])),4);
end