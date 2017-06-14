function yt =  operatorTVTranspose(y, T)

yt=real(ifft2(conj(T(:,:,1)).*fft2(y(:,:,1)))) + real(ifft2(conj(T(:,:,2)).*fft2(y(:,:,2))));

    