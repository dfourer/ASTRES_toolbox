function e = eulerian(k, n)
% function e = eulerian(k, n)
% 
% Compute the Eulerian number
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

e = 0;

for j = 0:n
  e = e + (-1)^j * binomial(k+1, j) * (n+1-j)^k;
end

end
