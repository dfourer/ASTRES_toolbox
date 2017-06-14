function [tfr, ma] = tfrstzp(x, M, w0T, gamma_K)
% [tfr] = tfrstzp(x, M, w0T, gamma_K)
% Compute the discrete-time S-transform and store energy before and after
% the transform for exact reconstruction using rectfrst()
% Computation time is significantly higher.
% 
% Usage example:
% [tfr, before, after] = tfrstzp(x, M, w0T, gamma_K);
% [x_hat] = rectfrst(tfr)
% x_hat   = x_hat(2:(end-1))    %% ignore before and after sample
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

before = zeros(M, 1);
after  = zeros(M, 1);

x_ref = x;
%% apply zero-padding
z_pad  = round(w0T * M * sqrt(2*log(1/gamma_K)) / pi);   %% maximal width
x = [zeros(1, z_pad)  x_ref  zeros(1, z_pad)];
i_s = (z_pad+1):(length(x)-z_pad);

tfr = zeros(M, N);

K = w0T * M * sqrt(2*log(1/gamma_K)) / pi;
L = sqrt(2*pi) / (M  * w0T);
A = -2 * pi^2 / ((M*w0T)^2);
B = -1i * 2*pi / M;

ma = m_axis(M); % m axis
for n = 1:length(x) %%N
  for m = 1:M
    %nn = n-1;
    
    nn = n-i_s(1)  + 1; %% dephase transform signal of 1 for reconstruction
    
    mm = ma(m);
    m2 = mm^2;
    
    K_m = K/abs(mm);
    L_m = abs(mm) * L;
    A_m2 = A * m2;
    B_m  = B * mm;
    
    k_min = min(n-1, round(K_m/2));
    k_max = min(length(x)-n, round(K_m/2));

    k = (-k_min):k_max;
    k2 = k.^2;
    
    tfr_tmp = L_m * exp(B_m * nn) * sum( x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .* k));
    %tfr(n,m) = L_m * exp(B_m * nn) * sum( x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .* k));
    
    if n < i_s(1)
      before(m) = before(m) + tfr_tmp;
      %fprintf(1, '1\n')
    elseif n > i_s(end)
      after(m) = after(m) + tfr_tmp;
      %fprintf(1, '2\n')
    else
      %fprintf(1, '3\n')
      tfr(m, n-i_s(1)+1)    = tfr_tmp;
    end
    
  end
end

%% add marginalized over time transform before and after the tfr
tfr  = [before tfr after];
end