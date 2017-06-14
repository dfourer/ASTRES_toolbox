function value = SNR(s, s_hat)
% function value = SNR(s, s_hat)
%
% Compute the signal to noise ratio expressed in dB
% defined by 10 log10(|s|^2 / |s-s_hat|^2)
%
%  INPUT:  
%  s:       Reference signal
%  s_hat:   Estimated signal
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

value=20*log10(norm(s)/norm(s-s_hat));

end