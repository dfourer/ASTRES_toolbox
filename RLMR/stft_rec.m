function [x_hat] = stft_rec(tfr, k, L, mi, mf, M, n0)
% [x_hat] = stft_rec(tfr, k, L, m, M, n0)
%
% Signal reconstruction from recursive STFT
%
% INPUT:
% tfr:  recursive STFT
% k:    order
% L:    window length
% mi, mf: mi,mf:        initial and final frequency bin in range [0,M]
% M:    number of frequency bin
% n0:   time shift
%
% OUTPUT:
% x_hat: reconstructed signal
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Feb. 2016
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

default_n0 = (k-1) * L;
m          = mi:mf;


if ~exist('n0', 'var')
  n0 = default_n0;
end

if length(m) ~= size(tfr, 1)
 error('Invalid values for mi and mf')    
end

if sum(sum(abs(tfr(:,1:n0)))) > eps
 warning('Energy lost during the reconstruction, try a lower value for n0');  
end

h02 = n0^(k-1)/(L^k*factorial(k-1)) * exp(-n0/L);   %% h0 * Ts
delta_w = 1 / M;                                    %% delta_w * Ts

if h02 == 0
  warning(1, 'Invalid value for n0, use default (k-1)*L ');
  n0 = default_n0;
end

%% (1)
% tmp = zeros(1,size(tfr, 2));
% for n = 1:size(tfr, 2)
%  for i = 1:length(m)
%   tmp(n) = tmp(n) + tfr(i,n) * exp(-2*pi*1i*m(i)/M * n0);   
%  end
% end

%% faster than (1)
x_hat = delta_w/h02 * sum(tfr .* repmat(exp(-2*pi*1i* m'/M * n0),1, size(tfr,2)), 1);

%% shift reconstructed signal
x_hat = x_hat((n0+1):end);
x_hat = x_hat(:);
end