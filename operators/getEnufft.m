function [E,out] = getEnufft(sense_maps,Nt,varargin)
% uniform sampling case: considers basis, fmap, nufft
% cl: considers nufft based on Otazo's case
arg.Nx = size(sense_maps,1);
arg.Ny = size(sense_maps,2);
arg.Nt = Nt;
arg.Nc = size(sense_maps,3);
arg.mask = true(arg.Nx, arg.Ny, arg.Nt);
arg.smaps = sense_maps;
arg.samp = []; 
arg.N = [arg.Nx,arg.Ny];
arg.fov = [22,22];
arg.basis = {'dirac'};
arg.ksp = mri_trajectory('cartesian', {}, arg.N, arg.fov);
arg.M = size(arg.ksp,1);
%nufft
arg.om = [];
arg.wi = [];%cl: this doesn't come up in simulation, but otazo used it
arg.donufft = 0;
arg.Jd = [6,6];
arg.Kd = floor([arg.N*1.5]);
arg.n_shift = arg.N/2;
arg = vararg_pair(arg, varargin);
 
out.basistransform = E_basis(arg);
if arg.donufft
    arg.basistransform = out.basistransform;
else
    arg.basistransform = reshape(out.basistransform,[arg.N,arg.Nt]);
end

%nufft
if arg.donufft
    %input ksp, om w/ size [M,d,nt], and wi w/ size[M,nt]
    if isempty(arg.om)
        error('cl: need ksp,om,wi for nufft')
    else %construct st for nufft
        if size(arg.ksp,3) ~= arg.Nt
            error('cl: double check ksp dimension')
        end
        for tt=1:size(arg.ksp,3)
            arg.st{tt} = nufft_init(arg.om(:,:,tt), arg.N, arg.Jd, arg.Kd, arg.n_shift,'kaiser');
        end
    end
    %nufft size
    if (arg.Nc > 1)
        E = fatrix2('arg',arg,'imask', arg.mask, ...
                             'odim', [arg.M arg.Nt arg.Nc], 'forw', ...
                            @E_forw, 'back', @E_back);
    else 
        E = fatrix2('arg',arg ,'imask', arg.mask,...
                             'odim', [arg.M arg.Nt], 'forw', ...
                            @E_forw, 'back', @E_back);
    end
else
    %cartesian size
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
end

function transform = E_basis(arg)
switch arg.basis{1}
case 'dirac'
	Bi = ones(arg.M, arg.Nt);
case {'sinc', 'dirac*dx'}
	Bi = ones(arg.M, arg.Nt) * prod(arg.fov ./ arg.N);
case 'rect' % rect(x/dx) <=> |dx| * sinc(dx * u)
	dx = abs(arg.fov) ./ arg.N; % usual default
	Bi = ones(arg.M,arg.Nt);
	for id=1:size(arg.ksp,2)
		if dx(id)
			Bi = Bi .* (nufft_sinc(dx(id) * squeeze(arg.ksp(:,id,:))));
		end
    end
otherwise
	fail('unknown basis_type "%s"', arg.basis)
end
transform = Bi; 
end


% y = G * x
function S = E_forw(arg,x)
if arg.donufft 
    S = zeros(arg.M,arg.Nt,arg.Nc);
    for tt=1:arg.Nt %parfor
        tmp = bsxfun(@times,x(:,:,tt),arg.smaps);
        S(:,tt,:) = reshape(nufft(tmp,arg.st{tt})/sqrt(prod(arg.N)),[arg.M,1,arg.Nc]);
    end
else 
    S = bsxfun(@times,x,reshape(arg.smaps,[arg.Nx,arg.Ny,1,arg.Nc]));
    S=fft2c_mri(S);
    if ~isempty(arg.samp) %cl: samp mask only when cartesian samp
        S = bsxfun(@times,S,arg.samp);
    end
end
S = bsxfun(@times,S,arg.basistransform);
end

% x = G' * y
function x = E_back(arg,S)
x = zeros(arg.Nx,arg.Ny,arg.Nt);
S = bsxfun(@times,S,conj(arg.basistransform));
if arg.donufft 
    if ~isempty(arg.wi)
        S = bsxfun(@times,S,arg.wi); %cl: from otazo
    end
    for tt=1:arg.Nt %cl: '/sqrt(prod(a.imSize))' from otazo
        tmp = reshape(nufft_adj(squeeze(S(:,tt,:)),arg.st{tt})/sqrt(prod(arg.N)),[arg.Nx,arg.Ny,arg.Nc]);
        x(:,:,tt)=sum(tmp.*conj(arg.smaps),3);
    end
else 
    if ~isempty(arg.samp) 
        S = bsxfun(@times,S,arg.samp);
    end
    s = ifft2c_mri(S);
    x = sum(bsxfun(@times,s,reshape(conj(arg.smaps),[arg.Nx,arg.Ny,1,arg.Nc])),4);
end
    
end