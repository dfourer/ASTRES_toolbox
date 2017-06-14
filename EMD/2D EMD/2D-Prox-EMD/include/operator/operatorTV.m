function yt =  operatorTV(y, T)

[n,m] = size(y);
yt = zeros(n,m,2);
yt(:,:,1)=real(ifft2(T(:,:,1).*fft2(y)));
yt(:,:,2)=real(ifft2(T(:,:,2).*fft2(y)));


    