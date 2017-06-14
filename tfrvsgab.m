function [tfr, stfr, lost, q_hatmap, if_hatmap] = my_tfrvsgab(x, M, L, q_method, if_method, gamma_K, q_threshold)
% [tfr] = my_tfrgab(x, M, L, gamma_K)
% Compute the discrete-time Levenberg-Marquard synchrosqueezed Gabor Transform
% 
% INPUT:
% x          : the signal to process
% M          : number of frequency bins to process (default: length(x))
% L          : window duration parameter:  w0 * T, (default: 10)
% q_method   : choose the chirp rate estimation method (default: 1 CRE1 , 2: CRE2t, 3: CRE2w, 4: CRE2r)
% if_method  : choose the IF estimation method (default: 1 biased , 2: unbiased, 3: reassignment operator)
% gamma_K    : threshold applied on window values (default: 10^(-4))
% q_threshold: threshold applied on window values (default: 10^(-4))
%
%
% OUTPUT:
% tfr       : STFT
% stfr      : synchrosqueezed STFT
% lost      : lost energy
% q_hatmap  : continuous time estimated CRE (\hat q)
% if_hatmap : continuous time estimated IF (\hat\omega)
%
% Author: D.Fourer
% Date: 27-06-2016
% Ref: [D. Fourer, F. Auger, K.Czarnecki and S. Meignen, Chirp rate and instantaneous frequency estimation. Proc. IEEE ICASSP 2017]

x = x(:).';          %convert x as a row vector
N = length(x);

if ~exist('M', 'var')
 M = N;
end
if ~exist('L', 'var')
 L = 10;
end
if ~exist('q_method', 'var')
 q_method = 2;
end
if ~exist('if_method', 'var')
 if_method = 1;
end
if ~exist('gamma_K', 'var')
 gamma_K = 10^(-4);
end

if ~exist('q_threshold', 'var')
 %% CRE threshold
 q_threshold = 1e-4; %eps
end


lost = 0;
tfr      = zeros(M, N);
stfr     = zeros(M, N);
q_hatmap = zeros(M, N);
if_hatmap= zeros(M, N);

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
  tdg_Ts   = -k.^2/(L^2) .* g_Ts;            %-t2g_Tsm1/(L^2);
  d2g_Ts3  = (-1/(L^2) + k.^2/L^4) .* g_Ts;  %-g_Ts/(L^2) +t2g_Tsm1/L^4; %(-1/(L^2) + k.^2/L^4) .* g_Ts;
  
  %% for q_method == 4
  t2dg     = k.^2 .* dg_Ts2;
  t3g_Tsm2 = -k .* t2g_Tsm1;
  
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
     
     if q_method == 4
      tfr_t2d_Tsm1 = exp(B * mm(m) * nn) * sum( x(n+k) .* t2dg     .* exp(B .* mm(m) .* k));
      tfr_t3_Tsm3  = exp(B * mm(m) * nn) * sum( x(n+k) .* t3g_Tsm2 .* exp(B .* mm(m) .* k));
      %tfr_dt       = tfr(m,n) + tfr_td;
      tfr_dt2_Tsm1 = 2 * tfr_t_Tsm1 + tfr_t2d_Tsm1;
      %tfr2         = tfr(m,n)^2;
     end
     
%     n_hat = n-round(real(tfr_t_Tsm1 / tfr(m,n)));
%     m_hat = m+round(M/(2*pi)* imag(tfr_d_Ts   /  tfr(m,n)) );  
     
%      %% reassigned coordinates
      n_tilde = n - round(tfr_t_Tsm1 / tfr(m,n));
      m_tilde = 1j * m + round(M/(2*pi) * tfr_d_Ts   /  tfr(m,n));
      
      n_tilde2 = n - tfr_t_Tsm1 / tfr(m,n);
      m_tilde2 = 1j * m + M/(2*pi) * tfr_d_Ts   /  tfr(m,n);
      
      if if_method == 1
        n_hat = round(real(n_tilde));
        m_hat = round(imag(m_tilde));
        n_hat2  = real(n_tilde2);
        m_hat2  = imag(m_tilde2);
      end
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRE1
     if q_method == 1 
       alpha_hat_denum = imag(tfr_t_Tsm1 * conj(tfr(m,n)));
       
       if abs(alpha_hat_denum) > q_threshold
         alpha_hat = real(tfr_d_Ts * conj(tfr(m,n))) / alpha_hat_denum;
       else
         alpha_hat = 0;
       end
       
       q_hat = 1j * alpha_hat;
       
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRE2t
     elseif q_method == 2 %% 

       q_hat_denum = tfr_t_Tsm1 * tfr_d_Ts - tfr_td * tfr(m,n);
       if abs(q_hat_denum) > q_threshold
         q_hat = (tfr_d2_Ts2 * tfr(m,n) - tfr_d_Ts^2) / q_hat_denum;
         
         %% d/dt(m_hat) / d/dt(n_hat)
         %q_hat = -imag( ST_D2g_Ts2_nm / ST_nm - (ST_Dg_Ts_nm/ST_nm)^2) / (1+ real( ST_TDg_nm / ST_nm - (ST_Tg_Tsm1_nm * ST_Dg_Ts_nm)/ST_nm^2));
       else
         q_hat = 0;
       end
       alpha_hat = imag(q_hat);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRE2w
     elseif q_method == 3
       q_hat_denum = tfr_t_Tsm1^2 - tfr_t2_Tsm2 * tfr(m,n);
       if abs(q_hat_denum) > q_threshold
         q_hat = ( tfr_td *tfr(m,n) -  tfr_t_Tsm1 * tfr_d_Ts + tfr(m,n)^2) / q_hat_denum;
       else
         q_hat = 0;
       end
       
       alpha_hat = imag(q_hat);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CRE2r
     elseif q_method == 4
       %a1 = (tfr_t_Tsm1 * tfr_d_Ts - tfr_td * tfr(m,n)) / tfr2;            % 1 - ((tfr_dt * tfr(m,n) - tfr_t_Tsm1 * tfr_d_Ts) / tfr2);  %  1 - d/dt[F_x^Th / F_x^h]
       %a2 = 1j * (tfr_t_Tsm1^2 - tfr_t2_Tsm2 * tfr(m,n)) / tfr2;
       %b1 = (tfr_t2d_Tsm1 * tfr(m,n) - tfr_t2_Tsm2 * tfr_d_Ts) / tfr2;     % (tfr_dt2_Tsm1 *  tfr(m,n) - tfr_t2_Tsm2 * tfr_d_Ts - 2 * tfr_t_Tsm1 *tfr(m,n))/ tfr2;
       %b2 = 1j * (tfr_t3_Tsm3 * tfr(m,n) - tfr_t2_Tsm2 * tfr_t_Tsm1) / tfr2;
       %y1 = (tfr_d2_Ts2 * tfr(m,n) - tfr_d_Ts^2) / tfr2;
       % y2 = 1j * (tfr_td *tfr(m,n) - tfr_d_Ts * tfr_t_Tsm1 - tfr2) / tfr2;
       

%        a1 = (tfr_t_Tsm1 * tfr_d_Ts - tfr_td * tfr(m,n));% / tfr2;            % 1 - ((tfr_dt * tfr(m,n) - tfr_t_Tsm1 * tfr_d_Ts) / tfr2);  %  1 - d/dt[F_x^Th / F_x^h]
%        a2 = (tfr_t_Tsm1^2 - tfr_t2_Tsm2 * tfr(m,n));% / tfr2;
%        b1 = (tfr_t2d_Tsm1 * tfr(m,n) - tfr_t2_Tsm2 * tfr_d_Ts);% / tfr2;     % (tfr_dt2_Tsm1 *  tfr(m,n) - tfr_t2_Tsm2 * tfr_d_Ts - 2 * tfr_t_Tsm1 *tfr(m,n))/ tfr2;
%        b2 = (tfr_t3_Tsm3 * tfr(m,n) - tfr_t2_Tsm2 * tfr_t_Tsm1);% / tfr2;
%        y1 = (tfr_d2_Ts2 * tfr(m,n) - tfr_d_Ts^2);% / tfr2;
%        y2 = (tfr_td *tfr(m,n) - tfr_d_Ts * tfr_t_Tsm1 - tfr2);% / tfr2;
       
       
       Ar = [tfr_t2_Tsm2  -tfr_t_Tsm1  tfr(m,n);...
            tfr_dt2_Tsm1  -tfr_td      tfr_d_Ts;...
            tfr_t3_Tsm3  -tfr_t2_Tsm2  tfr_t_Tsm1];      
      %tfr_t2d_Tsm1
      % b2 = 1j * (Tsm3_yT3g * tfr_k(m,n) - Tsm2_yT2g * Tsm1_yTg) / tfr_k(m,n)^2;
      
      % b1 = (Tsm1_yDT2g * tfr_k(m,n) - Tsm2_yT2g * Ts_yDg) / tfr_k(m,n)^2;

      % y1 = (Ts2_yD2g * tfr_k(m,n) - Ts_yDg^2) / tfr_k(m,n)^2;
      % y2 = 1j * ((yTDg *  tfr_k(m,n)) - Ts_yDg * Tsm1_yTg) / tfr_k(m,n)^2;  
       
       %q_hat_denum = b2 * a1 - b1* a2;
       q_hat_denum = det(A);
       
       if abs(q_hat_denum) > q_threshold
         ur = [tfr_d_Ts tfr_d2_Ts2 tfr_td+tfr(m,n)].';
         Ar_inv = pinv(Ar);
         rx = Ar_inv(1,:) * ur;
         q_hat = Ar_inv(2,:) * ur - 2 * rx * (n-n_tilde);
         
         if abs(q_hat) > 1/q_threshold
           q_hat = 0;
         end
       else
         q_hat = 0;
       end
       
       alpha_hat = imag(q_hat);
     end

     %dphi_dw_Tsm1   = n-real( tfr_t_Tsm1  / tfr(m,n));
     %dphi_dt_Ts     = imag( tfr_d_Ts    /  tfr(m,n));
     %d2phi_dtdw     = real(tfr_td/tfr(m,n)       - ((tfr_t_Tsm1 * tfr_d_Ts)/tfr(m,n))^2  );
     %d2phi_dt2_Ts2  = imag(tfr_d2_Ts2/tfr(m,n)   - (tfr_d_Ts/tfr(m,n))^2);
     %d2phi_dw2_Tsm2 = -imag(tfr_t2_Tsm2/tfr(m,n) - (tfr_t_Tsm1/tfr(m,n))^2);
     q_hatmap(m, n) = q_hat;
     
     %% IF 1, biased when Ax is not constant
     if if_method == 1
       m_hat_q        = m_hat + round(M / (2*pi) * alpha_hat * (n-n_hat));
       if_hatmap(m,n) = ((m_hat2 + M / (2*pi) * alpha_hat * (n-n_hat2))-1)/M;
     elseif if_method == 2
     %% IF 2, unbiased
       m_hat_q        = imag( m_tilde + round(M / (2*pi) * q_hat * (n-n_tilde))); 
       if_hatmap(m,n) = (imag((m_tilde2 + M / (2*pi) * q_hat * (n-n_tilde2)))-1)/M;
     elseif if_method == 3
      %% IF 3, (classical IF (reassignment operator))
       m_hat_q        = imag(m_tilde);
       if_hatmap(m,n) = (imag(m_tilde2)-1)/M;
     end
     
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