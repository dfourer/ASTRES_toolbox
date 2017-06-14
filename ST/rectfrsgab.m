function [ s_hat ] = rectfrsgab(stfr, L, M)
% [ s_hat ] = recsMW( stfr, w0. T)
%
% signal reconstruction from the synchrosqueezed Gabor Transform
% formala
%
% INPUT:
% stfr     : the synchrosqueezed Gabor Transform
% L        : T * Ts, the length of the Gaussian Window
% M        : number of frequency bins
%
% OUTPUT:
% s_hat  : reconstructed signal
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]
% Ref: [P. Flandrin. Time-frequency/Time-scale analysis.  Wavelet analysis and its applications. Academic Press. vol 10. 1998]

s_hat = sum(stfr) * (2*pi)^(3/2) * L/M;

end

