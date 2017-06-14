function [ v ] = dhk( t, k, T)
% function [ v ] = dhk( t, k, T)
%
% Return the value of the first derivative of the special case window
%
% INPUT:
% t:    time axis
% k:    order
% T:    window length
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if k < 1
  error('k should be greater or equal to 1')
elseif k == 1
 v = -T^(-2) * exp(-t/T);
else 
 v = T^(-1) * (hk(t,k-1,T) - hk(t, k, T));
end

end

