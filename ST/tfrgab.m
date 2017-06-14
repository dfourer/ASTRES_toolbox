function [tfr] = tfrgab(x, M, L, gamma_K)
% [tfr] = my_tfrgab(x, M, L, gamma_K)
% Compute the discrete-time Gabor Transform
% 
% INPUT:
% x      : the signal to process
% M      : number of frequency bins to process (default: length(x))
% L      : window duration parameter:  w0 * T, (default: 10)
% gamma_K: threshold applied on window values (default: 10^(-4))
%
% OUTPUT:
% tfr    : discrete Stockwellogram
%
% Author: D.Fourer
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. Signal Processing. 2015]

x = x(:).';          %convert x as a row vector
N = length(x);

if ~exist('M', 'var')
 M = N;
end
if ~exist('L', 'var')
 L = 10;
end
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end

tfr = zeros(M, N);


K = 2 * L * sqrt(2*log(1/gamma_K));  %% window length in samples


A = 1/(sqrt(2*pi)*L);
B = -1i * 2*pi / M;
C = -1 / (2*L^2);

mm = m_axis(M);
for n = 1:N
  
  k_min = min(n-1, round(K/2));
  k_max = min(N-n, round(K/2));
  k = (-k_min):k_max;
  k2 = k.^2;
  g = A * exp( C * k2);

%   z = zeros(1, M-length(g));
%   mm = m_axis(M);
%   tfr(n,:) = fft([z x(n+k) .* g], M) .* exp(B * mm * (n-1));

  for m = 1:M
    tfr(m,n) = exp(B * mm(m) * (n-1)) * sum( x(n+k) .* g .* exp(B .* mm(m) .* k));
  end
end

%% used to obtain matlab conventions
%tfr  = transpose(tfr);
end