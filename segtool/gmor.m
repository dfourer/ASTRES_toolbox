function psih = gmor(beta,gamma,xi)
% gmor : the generalized Morse wavelet in Fourier domain
    psih = 2*(exp(1)*gamma/beta)^(beta/gamma) * xi.^beta .* exp(-xi.^gamma);
    psih(xi<0)=0;
end