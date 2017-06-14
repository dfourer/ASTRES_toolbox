function [tfr, ma] = tfrst(x, M, w0T, gamma_K)
% [tfr] = tfrst(x, M, w0T, gamma_K)
% Compute the discrete-time S-transform
% 
% INPUT:
% x      : the signal to process
% M      : number of frequency bins to process (default: length(x))
% w0T    : central frequency * window duration parameter:  w0 * T, (default: 2pi 50)
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
if ~exist('w0T', 'var')
 w0T = 2*pi*50;
end
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end

tfr = zeros(M, N);

K = w0T * M * sqrt(2*log(1/gamma_K)) / pi;
L = sqrt(2*pi) / (M  * w0T);
A = -2 * pi^2 / ((M*w0T)^2);
B = -1i * 2*pi / M;

ma = m_axis(M); % m axis
%Mh = floor(M/2);
for n = 1:N
  for m = 1:M
    nn = n-1;
%     if m <= Mh   %% map to [0, M/2]
%      mm = m-1;
%     else         %% map to [-M/2, 0]
%      mm = -M+m;
%     end
    mm = ma(m);
    m2 = mm^2;
    
    K_m = K/abs(mm);
    L_m = abs(mm) * L;
    A_m2 = A * m2;
    B_m  = B * mm;
    
    k_min = min(nn , round(K_m/2));
    k_max = min(N-n, round(K_m/2));

    k = (-k_min):k_max;
    k2 = k.^2;
    tfr(m,n)    = L_m * exp(B_m * nn) * sum( x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .* k));
  end
end

end