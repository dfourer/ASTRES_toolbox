function [ s_hat ] = recsMW( stfr, w0, T, m_range)
% [ s_hat ] = recsMW( stfr, w0. T)
%
% signal reconstruction from synchrosqueezed Morlet wavelet transform
%
% INPUT:
% stfr     : synchrosqueezed Morlet wavelet transform
% w0       : central frequency in radian of the Morlet wavelet transform
% T        : T parameter of the Morlet wavelet transform (related to the length of a Gaussian window)
% m_range  : frequency bins to integrate during the reconstruction (default 1:M)
%
% OUTPUT:
% s_hat  : reconstructed signal
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]
% Ref: [I. Daubechies, J. Lu and H-T Wu, Synchrosqueezed wavelet transforms: An empirical mode decomposition-like tool,
%       Applied and Computational Harmonic Analysis. vol 30(2), pages 243-261. 2011]

if(~exist('m_range', 'var'))
 m_range = 1:size(stfr, 1);
end

R = sum(stfr(m_range,:));


K = 500;
dw = 0.001;
w = (w0-K):dw:(w0+K)+eps;

C_psi = sqrt(2*T * sqrt(pi))  * sum( 1./w .* exp(-((w0-w).^2 * T^2)/2)) * dw;

s_hat = 2/C_psi * R';

end

