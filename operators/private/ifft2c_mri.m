function x=ifft2c_mri(X)
% Claire Lin: Otazo's ifft
x=fftshift(fft(ifftshift(X,1),[],1),1)/sqrt(size(X,1));
x=fftshift(fft(ifftshift(x,2),[],2),2)/sqrt(size(X,2));
