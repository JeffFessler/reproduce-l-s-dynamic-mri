function T = getT(Nx,Ny,Nt,varargin)
arg.Nx = Nx;
arg.Ny = Ny;
arg.Nt = Nt;
arg.dim = 3;
arg = vararg_pair(arg, varargin);

T = fatrix2('idim', [arg.Nx arg.Ny arg.Nt], 'arg', ...
                        arg, 'odim', [arg.Nx arg.Ny arg.Nt], 'forw', ...
                        @T_forw, 'back', @T_back);
end


% y = G * x
function S = T_forw(arg,b)
S = fftshift(fft(b,[],arg.dim),arg.dim)/sqrt(size(b,arg.dim)); 
end

% x = G' * y
function x = T_back(arg,b)
x = ifft(ifftshift(b,arg.dim),[],arg.dim)*sqrt(size(b,arg.dim));
end