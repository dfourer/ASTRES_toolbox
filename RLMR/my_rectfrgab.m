function [ s_hat ] = my_rectfrgab(tfr, L, M, mm)
% [ s_hat ] = my_recsMW( stfr, w0. T)
%
% signal reconstruction from the Gabor Transform using the simplified
% formula
%
% INPUT:
% tfr      : the synchrosqueezed Gabor Transform
% L        : T * Ts, the length of the Gaussian Window
% M        : number of frequency bins
% mm       : frequencies bins to reconstruct
%
% OUTPUT:
% s_hat    : reconstructed signal
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]
% Ref: [P. Flandrin. Time-frequency/Time-scale analysis.  Wavelet analysis and its applications. Academic Press. vol 10. 1998]

N = size(tfr, 2);

if ~exist('mm', 'var')
 mm = m_axis(M);
end

s_hat = zeros(1, N);
for n = 1:N
 
 s_hat(n) = sum(tfr(:, n) .* exp(1i *2*pi*(n-1)*mm'/M));
end

s_hat = s_hat * sqrt(2*pi) * L/M;

end

