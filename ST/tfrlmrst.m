function [tfr, rtfr, lost] = tfrlmrst(x, mu, M, w0T, gamma_K)
% [tfr, rtfr, lost] = tfrrst(x, M, w0T, gamma_K)
% Compute the discrete Levenberg-Maquardt Reassigned S transform
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
% lost   : amount of lost energy during reassignment
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

tfr     = zeros(N, M);
rtfr    = zeros(N, M);

tfr_t_Tsm1   = zeros(N, M);    %% ST_x^{Th}   * Ts^{-1}
tfr_t2_Tsm2  = zeros(N, M);    %% ST_x^{T^2h} * Ts^{-2}
tfr_t3_Tsm3  = zeros(N, M);    %% ST_x^{T^3h} * Ts^{-3}


K = w0T * M * sqrt(2*log(1/gamma_K)) / pi;
L = sqrt(2*pi) / (M * w0T);
A = -2 * pi^2 / ((w0T)^2*M^2);
B = -1i * 2*pi / M;

Ld = (2*pi)^2 / (w0T^2 * M^2);
%Mh = round(M/2);
ma = m_axis(M); % m axis
%% Compute the classical Stockwellogram
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
    
    K_m  = K/abs(mm);
    L_m  = abs(mm) * L;
    A_m2 = A * m2;
    B_m  = B * mm;
    
    k_min = min(n-1, round(K_m/2));
    k_max = min(N-n, round(K_m/2));
    k = (-k_min):k_max;
    k2 = k.^2;
    
    C_mn = L_m * exp(B_m * nn);
    D_mk = x(n+k) .* exp( A_m2 .* k2) .* exp(B_m .*k);
    tfr(n,m)           =  C_mn * sum(         D_mk);
    tfr_t_Tsm1(n,m)    = -C_mn * sum( k    .* D_mk);
    tfr_t2_Tsm2(n,m)   =  C_mn * sum( k.^2 .* D_mk);
    tfr_t3_Tsm3(n,m)   = -C_mn * sum( k.^3 .* D_mk);

    %% Compute the reassigned Stockwellogram
    if abs(tfr(n,m)) > eps
      
      %% 1 computations of ST_x
      ST_nm          = tfr(n,m);
      ST_Tg_Tsm1_nm  = tfr_t_Tsm1(n,m);
      ST_T2g_Tsm2_nm = tfr_t2_Tsm2(n,m);
      ST_T3g_Tsm3_nm = tfr_t3_Tsm3(n,m);
      
      ST_Dg_Ts_nm    = -Ld * m2 * ST_Tg_Tsm1_nm;
      ST_TDg_nm      = -Ld * m2 * ST_T2g_Tsm2_nm;
      ST_D2g_Ts2_nm  = -Ld * m2 * ( ST_nm + ST_TDg_nm); %Ld * m2 * ( ST_nm - ST_TDh_nm);
      
      omega_Tsm1 = 2*pi*abs(mm)/M;
      %% 2 computation of partial derivatives
      R      = [real(ST_Tg_Tsm1_nm / ST_nm); -imag(ST_Dg_Ts_nm / ST_nm)];
      
      dR1_dt     = 1 + real(ST_TDg_nm  / ST_nm - ((ST_Tg_Tsm1_nm * ST_Dg_Ts_nm) / ST_nm^2));
      dR2_dt_Ts2 = -imag(ST_D2g_Ts2_nm / ST_nm - (ST_Dg_Ts_nm   / ST_nm)^2 );
      dR1_dw_Tsm2= -imag(ST_T2g_Tsm2_nm / ST_nm - (ST_Tg_Tsm1_nm / ST_nm)^2 ) + omega_Tsm1/(w0T)^2 * real((ST_T2g_Tsm2_nm*ST_Tg_Tsm1_nm)/ST_nm^2-ST_T3g_Tsm3_nm/ST_nm);
      dR2_dw     = -real(ST_TDg_nm      / ST_nm - ((ST_Tg_Tsm1_nm * ST_Dg_Ts_nm) / ST_nm^2)) - omega_Tsm1^3/(w0T)^4 * imag(ST_T3g_Tsm3_nm/ST_nm - (ST_T2g_Tsm2_nm*ST_Tg_Tsm1_nm)/ST_nm^2)...
                   +2*omega_Tsm1/(w0T)^2 * imag(ST_Tg_Tsm1_nm/ST_nm);
  
      %dR1_dt     = 1 + real(ST_TDg_nm  / ST_nm - ((ST_Tg_Tsm1_nm * ST_Dg_Ts_nm) / ST_nm^2));
      %dR2_dt_Ts2 = -imag(ST_D2g_Ts2_nm / ST_nm - (ST_Dg_Ts_nm   / ST_nm)^2 );
      %dR1_dw_Tsm2= -imag(ST_T2g_Tsm2_nm / ST_nm - (ST_Tg_Tsm1_nm / ST_nm)^2 );
      %dR2_dw     = -real(ST_TDg_nm      / ST_nm - ((ST_Tg_Tsm1_nm * ST_Dg_Ts_nm) / ST_nm^2));

      %NablaR = [dR1_dt dR1_dw_Tsm2; dR2_dt_Ts2 dR2_dw];R_hat = pinv(NablaR +  mu *  eye(2)) * R;
      %% about 2x faster
      NablaR_pmuI2 = [dR1_dt+mu dR1_dw_Tsm2;dR2_dt_Ts2 dR2_dw+mu];
      det_NR = NablaR_pmuI2(1,1)*NablaR_pmuI2(2,2)-NablaR_pmuI2(1,2)*NablaR_pmuI2(2,1);
      R_hat  = 1/det_NR * [NablaR_pmuI2(2,2) -NablaR_pmuI2(1,2); -NablaR_pmuI2(2,1) NablaR_pmuI2(1,1)] * R;
      
      
      n_hat = n - round(R_hat(1));
      m_hat = m - round(M / (2*pi) * R_hat(2));

      %% out of bounds
      n_out_of_bounds = n_hat < 1 || n_hat > N;
      m_out_of_bounds = m_hat < 1 || m_hat > M;
      
      if n_out_of_bounds || m_out_of_bounds
        lost = lost + abs(tfr(n,m))^2;
        continue;
      end
      
      rtfr(n_hat, m_hat) = rtfr(n_hat, m_hat) + abs(tfr(n,m))^2;
    end
  end
end

%% used to obtain matlab conventions
tfr  = transpose(tfr);
rtfr = transpose(rtfr);
end
