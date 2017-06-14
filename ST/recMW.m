function [ s_hat ] = recMW(Wx, w0, T, as)
% [ s_hat ] = recsMW( stfr, w0. T)
%
% signal reconstruction from Morlet wavelet transform using the truncated
% formala
%
% INPUT:
% Wx       : the Morlet wavelet transform of signal 
% w0       : central frequency in radian of the Morlet wavelet transform
% T        : T parameter of the Morlet wavelet transform (related to the length of a Gaussian window)
% as       : corresponding scales
%
% OUTPUT:
% s_hat  : reconstructed signal
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]
% Ref: [P. Flandrin. Time-frequency/Time-scale analysis.  Wavelet analysis and its applications. Academic Press. vol 10. 1998]

N = size(Wx,2);  %% signal_length

K = 500;
dw = 0.001;
w = (w0-K):dw:(w0+K)+eps;

C_psi = sqrt(2*T * sqrt(pi))  * sum( 1./w .* exp(-((w0-w).^2 * T^2)/2)) * dw;


ds = abs([as(2)-as(1) diff(as)]);

s_hat = real(1/C_psi * sum(Wx .* repmat(as',1, N).^(-3/2) .* repmat(ds',1, N)));


end

