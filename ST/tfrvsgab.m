function [tfr, stfr, lost] = my_tfrvsgab(x, M, L,gamma_K)
% [tfr] = my_tfrgab(x, M, L, gamma_K)
% Compute the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform
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


lost = 0;
tfr    = zeros(M, N);
stfr   = zeros(M, N);

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
  
  %t2g_Tsm1 = k.^2 .* g_Ts;
  tdg_Ts   = -k.^2/(L^2) .* g_Ts;            %-t2g_Tsm1/(L^2);
  d2g_Ts3  = (-1/(L^2) + k.^2/L^4) .* g_Ts;  %-g_Ts/(L^2) +t2g_Tsm1/L^4; %(-1/(L^2) + k.^2/L^4) .* g_Ts;
  
  % fft(x(n+k) .* g, M);
  for m = 1:M
    nn = n-1;
    
    tfr(m,n) = exp(B * mm(m) * nn) * sum( x(n+k) .* g_Ts .* exp(B .* mm(m) .* k));  % 
    
    if abs(tfr(m,n)) > eps
        
     tfr_t_Tsm1 = exp(B * mm(m) * nn) * sum( x(n+k) .* tg     .* exp(B .* mm(m) .* k));  % 
     tfr_d_Ts   = exp(B * mm(m) * nn) * sum( x(n+k) .* dg_Ts2 .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 
     
     tfr_td      = exp(B * mm(m) * nn) * sum( x(n+k) .* tdg_Ts   .* exp(B .* mm(m) .* k));  % 
     %tfr_t2_Tsm2 = exp(B * mm(m) * nn) * sum( x(n+k) .* t2g_Tsm1 .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 
     tfr_d2_Ts2  = exp(B * mm(m) * nn) * sum( x(n+k) .* d2g_Ts3  .* exp(B .* mm(m) .* k));  %exp(B * mm * nn) * 

     %dphi_dw_Tsm1   = n-real( tfr_t_Tsm1  / tfr(m,n));
     %dphi_dt_Ts     = imag( tfr_d_Ts    /  tfr(m,n));
     
     %d2phi_dtdw     = real(tfr_td/tfr(m,n)       - ((tfr_t_Tsm1 * tfr_d_Ts)/tfr(m,n))^2  );
     %d2phi_dt2_Ts2  = imag(tfr_d2_Ts2/tfr(m,n)   - (tfr_d_Ts/tfr(m,n))^2);
     %d2phi_dw2_Tsm2 = -imag(tfr_t2_Tsm2/tfr(m,n) - (tfr_t_Tsm1/tfr(m,n))^2);


     %% reassigned coordinates
     n_hat = n-round(real(tfr_t_Tsm1 / tfr(m,n)));
     m_hat = m+round(M/(2*pi)* imag(tfr_d_Ts   /  tfr(m,n)) );      
      
      
     %% frequency modulation     
     %dt_hat_dt = -((ST_TDg_nm * ST_nm) - (ST_Tg_Tsm1_nm * ST_Dg_Ts_nm));
     dt_hat_dt = tfr_td - tfr_t_Tsm1 * tfr_d_Ts;
     if abs(dt_hat_dt) > eps
       %% works better !
       %q_hat = real( (tfr_d2_Ts2 * tfr(m,n) - tfr_d_Ts^2) / ( tfr(m,n)^2 +  (tfr_t_Tsm1 * tfr_d_Ts) - (tfr_td *tfr(m,n)) ));
              
       q_hat = imag( (tfr_d2_Ts2 * tfr(m,n) - tfr_d_Ts^2) / (  (tfr_t_Tsm1 * tfr_d_Ts) - (tfr_td *tfr(m,n)) ));
       
       %% d/dt(m_hat) / d/dt(n_hat)
       %q_hat = -imag( ST_D2g_Ts2_nm / ST_nm - (ST_Dg_Ts_nm/ST_nm)^2) / (1+ real( ST_TDg_nm / ST_nm - (ST_Tg_Tsm1_nm * ST_Dg_Ts_nm)/ST_nm^2));
      else
       q_hat = 0;
     end

      %q_hat = 0;
      
      m_hat_q = m_hat + round(M / (2*pi) * q_hat * (n-n_hat));
     
     %% out of bounds
     m_out_of_bounds = m_hat_q < 1 || m_hat_q > M;
     
     if  m_out_of_bounds
      lost = lost + abs(tfr(m,n))^2;
      continue;
     end
      
     stfr(m_hat_q,n) = stfr(m_hat_q,n) + tfr(m,n)/(2*pi) * exp(2*1i*pi*mm(m)*nn/M);
    end
     
  end %% m
end %% n

%% used to obtain matlab conventions
%tfr  = transpose(tfr);
%stfr  = transpose(stfr);
end