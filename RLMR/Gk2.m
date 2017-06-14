function [a, b] = Gk2(k, L, alpha)
% [a, b] = Gk2(k, L, alpha)
%
% Compute coefficients of filter G_k(z,w)
%
% INPUT:
% k:        order
% L:        window length (samples)
% alpha:    exp(pTs)
%
% OUTPUT:
% filters coefficients
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

b = zeros(1, k);
a = zeros(1, k+1);

for i = 1:(k+1)
  if (i<=k)
   b(i) = eulerian(k-1, k-i) * alpha^(i-1)/factorial(k-1)/L^k;
  end  
  a(i) = (-alpha)^(i-1) * binomial(k, i-1);
end

end
