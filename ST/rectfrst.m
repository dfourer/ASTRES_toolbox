function [s_hat] = rectfrst(tfr, m, M)
% [s_hat] = rectfrst(tfr)
% 
% Signal reconstruction from S-transform (method 1)
% using classical inversion formula. 
% /!\ For better results, signal has to be zeropadded before computing
% transform as:
%
%%%%%%%%%%%%%%%% Usage Example %%%%%%%%%
%  %% 0-pad
%  z_pad   = round(w0T * M * sqrt(2*log(1/gamma_K)) / pi);
%  s0p     = [zeros(z_pad,1) ; s ; zeros(z_pad, 1)];
%  i_s     = z_pad+1:(length(s0p)-z_pad);
%
%  %% Compute transfirm
%  [tfr0p] =  tfrst(s0p, M, w0T, gamma_K);
%
%  %% reconstruction
%  stmp    = rectfrst(tfr0p);
%  s_hat   = stmp(i_s);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% tfr   : input S-transform
% m     : frequency bin to reconstruct in [-M/2, M/2]
% M     : number of frequency bins
%
% OUTPUT:
% s_hat : reconstructed signal
%
% Author: D.Fourer
% Date: 05-06-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. on Signal Processing. 2015]

N = size(tfr,2);

if ~exist('m', 'var')
  %warning('reconstruct with default frequency bins')
  M = size(tfr,1);
  m = m_axis(M);
end


% s_hat = zeros(1, N);
% 
% for n = 1:N
%  for m = 1:M
%   for k = 1:N
%    s_hat(n) = s_hat(n) + tfr(m, k) * exp(1i * 2 * pi * m * (n-1)/M) / M;
%   end
%  end
% end

F_x = transpose(sum(tfr,2));

%s_hat = ifft(F_x, N);
s_hat = zeros(1, N);
A = 2 * 1i * pi / M;
for n = 1:N
  s_hat(n) = 1/M * sum( F_x .* exp(A .* m .* (n-1)));
end


end