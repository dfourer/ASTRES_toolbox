function [x_hat] = sstft_rec(stfr, k, L, M, n0)
% [x_hat] = sstft_rec(stfr, k, L, M, n0)
%
% signal reconstruction from synchrosqueezed recursive STFT
%
% INPUT:
% stfr: the synchrosqueezed recursive STFT
% k: order
% L: window length
% M: number of frequency bin
% n0: time-shift used to allow reconstruction with a causal window (should be greater than 0, default:1)
%
% OUTPUT:
% x_hat: reconstructed signal
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if ~exist('n0', 'var')
 n0 = 1;    
end

if sum(sum(abs(stfr(:,1:n0)))) > eps
 warning('Energy lost during the reconstruction, try a lower value for n0')   
end

h02 = n0^(k-1)/(L^k*factorial(k-1)) * exp(-n0/L);   %% h0 * Ts
delta_w = 1 / M;                                    %% delta_w * Ts

x_hat = delta_w/h02 * sum(stfr, 1);

%% shift reconstructed signal
x_hat = x_hat((n0+1):end);
x_hat = x_hat(:);
end
