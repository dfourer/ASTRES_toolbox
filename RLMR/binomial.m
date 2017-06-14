function e = binomial( n, k )
% function e = binomial( n, k )
%
% Computes the binomial coefficient n!/(k!(n-k)!)
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

e = factorial(n) / (factorial(k) * factorial(n-k));

end

