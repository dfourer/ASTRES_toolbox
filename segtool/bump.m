function psih = bump(mu,sigma,xi)
% bump : the bump wavelet in Fourier domain
    psih = exp(1-1./(1-((xi-mu)./sigma).^2));
    psih(abs(xi-mu)>sigma)=0;
end