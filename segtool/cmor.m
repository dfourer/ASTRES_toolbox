function psih = cmor(Fb,Fc,xi)
% cmor : complex Morlet Wavelet in Fourier domain
    psih = exp(-Fb*pi^2*(xi-Fc).^2);
end