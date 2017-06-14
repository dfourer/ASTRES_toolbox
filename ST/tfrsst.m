function [tfr, stfr, lost] = tfrsst(x, M, w0T, gamma_K)
% [tfr, stfr, lost] = tfrsst(x, M, w0T, gamma_K)
% Compute the discrete-time synchrosqueezed S-transform
% 
% INPUT:
% x        : the signal to process
% M        : number of frequency bins to process (default: length(x))
% w0T      : central frequency * window duration parameter:  w0 * T, (default: 2pi 50)
% gamma_K  : threshold applied on window values (default: 10^(-4))
%
% OUTPUT:
% tfr      : discrete Stockwellogram
% stfr     : discrete Synchrosqueezed Stockwellogram
% lost     : amount of lost energy during synchrosqueezing
%
% Author: D.Fourer
% Date: 05-06-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. Signal Processing. 2015]

x = x(:).';          %convert x as a row vector
N = length(x);

lost = 0;

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
stfr    = zeros(M, N);
tfr_d   = zeros(M, N);


K = w0T * M * sqrt(2*log(1/gamma_K)) / pi;
L = sqrt(2*pi) / (M * w0T);
A = -2 * pi^2 / ((w0T)^2*M^2);
B = -1i * 2*pi / M;
D = -B; % 1i * 2* pi/ M;
ma = m_axis(M); % m axis
Ld = (2*pi)^2 / (w0T^2 * M^2);

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
    B_m  = B * mm;
    
    k_min = min(n-1, round(K_m/2));
    k_max = min(N-n, round(K_m/2));
    
    k = (-k_min):k_max;
    k2 = k.^2;
    tfr(m,n)      = L_m * exp(B_m * nn) * sum( x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .* k));
    tfr_d(m,n)    = Ld * m2 *L_m * exp(B_m * (n-1)) * sum( k .* x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .*k));
 
    %% Compute the synchrosqueezed Stockwellogram
    if abs(tfr(m,n)) > eps
      %n_hat = n - round(              real( tfr_t(m,n) / tfr(m,n) ) );
      m_hat = m + round( M / (2*pi) * imag( tfr_d(m,n) / tfr(m,n) ) );
      
      %% out of bounds (Should never occur)
      m_out_of_bounds = m_hat < 1 || m_hat > M;
      if m_out_of_bounds
        lost = lost + abs(tfr(m,n))^2;
        continue;
      end
      
      stfr(m_hat, n) = stfr(m_hat, n) +  1/abs(mm) * exp(D * (n-1) * mm) * tfr(m,n);
    end
  end
  stfr(:, n) = stfr(:, n) .* abs(ma).'; %((1:M)-1+eps);
end

end