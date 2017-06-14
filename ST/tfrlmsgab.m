function [tfr, rtfr, lost] = tfrlmsgab(x, M, L, mu, gamma_K)
% [tfr] = tfrgab(x, M, L, gamma_K)
% Compute the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform
% 
% INPUT:
% x      : the signal to process
% M      : number of frequency bins to process (default: length(x))
% L      : window duration parameter:  w0 * T, (default: 10)
% mu     : damping parameter 
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


lost = 0;
tfr    = zeros(M, N);
rtfr   = zeros(M, N);

% tfr_t = zeros(N, M);
% tfr_d = zeros(N, M);
% 
% tfr_td = zeros(N, M);
% tfr_d2 = zeros(N, M);
% tfr_t2 = zeros(N, M);


K = 2 * L * sqrt(2*log(1/gamma_K));  %% window length in samples

%T = Ts * L;

A = 1/(sqrt(2*pi)*L);
B = -1i * 2*pi / M;
C = -1 / (2*L^2);

%Mh = floor(M/2);
mm = m_axis(M);
for n = 1:N
  
  k_min = min(n-1, round(K/2));
  k_max = min(N-n, round(K/2));

  k = (-k_min):k_max;
  k2 = k.^2;
  g_Ts     = A * exp( C * k2);
  tg       = -k .* g_Ts;
  dg_Ts2   = -1/(L^2) .*tg;    %L^(-2) * k .* g_Ts;
  
  t2g_Tsm1 = k.^2 .* g_Ts;
  tdg_Ts   = -t2g_Tsm1/(L^2); %-k.^2/(L^2) .* g_Ts;            %-t2g_Tsm1/(L^2);
  d2g_Ts3  = -g_Ts/(L^2) + t2g_Tsm1/(L^4); %(-1/(L^2) + k.^2/L^4) .* g_Ts;  %-g_Ts/(L^2) +t2g_Tsm1/L^4; %(-1/(L^2) + k.^2/L^4) .* g_Ts;
  
  % fft(x(n+k) .* g, M);
  for m = 1:M
    nn = n-1;
    
    tfr(m,n) = exp(B * mm(m) * nn) * sum( x(n+k) .* g_Ts .* exp(B .* mm(m) .* k));  % 
    
    if abs(tfr(m,n)) > eps
        
     tfr_t_Tsm1 = exp(B * mm(m) * nn) * sum( x(n+k) .* tg     .* exp(B .* mm(m) .* k));  % 
     tfr_d_Ts   = exp(B * mm(m) * nn) * sum( x(n+k) .* dg_Ts2 .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 
     
     tfr_td      = exp(B * mm(m) * nn) * sum( x(n+k) .* tdg_Ts   .* exp(B .* mm(m) .* k));  % 
     tfr_t2_Tsm2 = exp(B * mm(m) * nn) * sum( x(n+k) .* t2g_Tsm1 .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 
     tfr_d2_Ts2  = exp(B * mm(m) * nn) * sum( x(n+k) .* d2g_Ts3  .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 

     dphi_dw_Tsm1   = -n+real( tfr_t_Tsm1  / tfr(m,n));
     dphi_dt_Ts     = imag( tfr_d_Ts    /  tfr(m,n));
     
     d2phi_dtdw     = real(tfr_td/tfr(m,n)       - ((tfr_t_Tsm1 * tfr_d_Ts)/tfr(m,n))^2  );
     d2phi_dt2_Ts2  = imag(tfr_d2_Ts2/tfr(m,n)   - (tfr_d_Ts/tfr(m,n))^2);
     d2phi_dw2_Tsm2 = -imag(tfr_t2_Tsm2/tfr(m,n) - (tfr_t_Tsm1/tfr(m,n))^2);
     
     R      = [n+dphi_dw_Tsm1; -dphi_dt_Ts];
     nablaR = [1+d2phi_dtdw+mu d2phi_dw2_Tsm2;-d2phi_dt2_Ts2 -d2phi_dtdw+mu];
     
     %faster than R_hat  = pinv(nablaR) * R;
     R_hat  = [nablaR(2,2) -nablaR(1,2);-nablaR(2,1) nablaR(1,1)] / (nablaR(1,1)*nablaR(2,2)-nablaR(2,1)*nablaR(1,2)) * R;
     
     
     %n_hat = n-round(R_hat(1)); % useless for the first order synchrosqueezing
     
     m_hat = m-round(M/(2*pi) * R_hat(2));
      
     %% out of bounds
     m_out_of_bounds = m_hat < 1 || m_hat > M;
     
     if  m_out_of_bounds
      lost = lost + abs(tfr(m,n))^2;
      continue;
     end
      
     rtfr(m_hat,n) = rtfr(m_hat,n) + tfr(m,n)/(2*pi) * exp(2*1i*pi*mm(m)*nn/M);
    end
     
  end %% m
end %% n

%% used to obtain matlab conventions
%tfr  = transpose(tfr);
%rtfr  = transpose(rtfr);
end