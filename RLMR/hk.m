function [ v ] = hk( t, k, T)
% [ v ] = hk( t, k, T)
%
% Return the value of the special case window
%
% INPUT:
% t: time axis
% k: filter-order
% T: duration in seconds
%
% OUTPUT:
% v:  values of h_k(t)
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if k < 1
  error('k should be greater or equal to 1')  
end

v = 1 / ((T^k) * factorial(k-1)) *  t.^(k-1) .* exp(-t/T);

end

