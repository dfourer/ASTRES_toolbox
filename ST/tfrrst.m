function [tfr, rtfr, ma, lost] = tfrrst(x, M, w0T, gamma_K)
% [tfr, rtfr, lost] = tfrrst(x, M, w0T, gamma_K)
% Compute the discrete-time Reassigned S transform
% 
% INPUT:
% x      : the signal to process
% M      : number of frequency bins to process (default: length(x))
% w0T    : central frequency * window duration parameter:  w0 * T, (default: 2pi 50)
% gamma_K: threshold applied on window values (default: 10^(-4))
%
% OUTPUT:
% tfr    : discrete Stockwellogram
% rtfr   : discrete Reassigned Stockwellogram
% ma     : frequency bins
% lost   : amount of lost energy during reassignment
%
% Author: D.Fourer
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. Signal Processing. 2015]

x = x(:).';          %convert x as a row vector
N = length(x);

lost = 0;            % lost energy after reassignment

if ~exist('M', 'var')
 M = N;
end
if ~exist('w0T', 'var')
 w0T = 2*pi*50;
end
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end

tfr     = zeros(M, N);
rtfr    = zeros(M, N);

tfr_t   = zeros(M, N);
tfr_d   = zeros(M, N);


K = w0T * M * sqrt(2*log(1/gamma_K)) / pi;
L = sqrt(2*pi) / (M * w0T);
A = -2 * pi^2 / ((w0T)^2*M^2);
B = -1i * 2*pi / M;

Ld = (2*pi)^2 / (w0T^2 * M^2);

ma = m_axis(M); % m axis
%Mh = floor(M/2);
for n = 1:N
  for m = 1:M
    nn = n-1;
%      if m <= Mh   %% map to [0, M/2]
%       mm = m-1;
%      else         %% map to [-M/2, 0]
%       mm = -M+m;
%      end
    mm = ma(m);
    m2 = mm^2;
    
    K_m = K/abs(mm);
    L_m = abs(mm) * L;
    A_m2 = A * m2;
    B_m = B * mm;
    
    k_min = min(nn , round(K_m/2));
    k_max = min(N-n, round(K_m/2));
    k = (-k_min):k_max;
    k2 = k.^2;
    
    tfr(m,n)      = L_m * exp(B_m * nn) * sum( x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .*k));
    
    tfr_t(m,n)    = -L_m * exp(B_m * nn) * sum( k .* x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .*k));
    tfr_d(m,n)    = -Ld * m2 * tfr_t(m,n);
    
    %% Compute the reassigned Stockwellogram
    if abs(tfr(m,n)) > eps
      n_hat = n - round(              real( tfr_t(m,n) / tfr(m,n) ) );
      m_hat = m + round( M / (2*pi) * imag( tfr_d(m,n) / tfr(m,n) ) );
      
      %% out of bounds
      n_out_of_bounds = n_hat < 1 || n_hat > N;
      m_out_of_bounds = m_hat < 1 || m_hat > M;
      
      if n_out_of_bounds || m_out_of_bounds
        lost = lost + abs(tfr(m,n))^2;
        continue;
      end
      
      rtfr(m_hat,n_hat) = rtfr(m_hat,n_hat) + abs(tfr(m,n))^2;
    end
  end
end

end
