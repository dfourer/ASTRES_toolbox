function [f, a_n,b_n,Freq] = get_fscoeff(y,nfft,t,SampFreq)
% To obtain parameters for GWarblet
% y : the extracted IF
% nfft: the number of fft

% f : reconstructed fourier series
% a_n : coeff of sin;
% b_n : coeff of cos;
yfft = fft(y,nfft)/nfft;
Nq = ceil(nfft/2);
Freq = SampFreq/nfft*[0:Nq-1];

% Build the fourier series
f = zeros(size(y));
for j = 2:round(Nq/3)
    a_n(j-1) = -2*imag(yfft(j));
    b_n(j-1) = 2*real(yfft(j));
    f = f + b_n(j-1)*cos(2*pi*t*Freq(j)) + a_n(j-1)*sin(2*pi*t*Freq(j));
end
f = abs(yfft(1))+f;