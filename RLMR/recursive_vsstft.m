function [tfr, stfr, nfreqs, below, over, q_hat_map, if_map] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_method, if_method, q_threshold)
% [tfr, stfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_method, if_method)
%
% Compute the recursive vertical synchrosqueezed STFT
%
% INPUT:
% x          : input signal
% k          : order
% L          : window length
% mi,mf      : initial and final frequency bin in range [0,M]
% M          : number of frequency bins
% n0         : time-shift used to allow reconstruction with a causal window (should be greater than 0, default: 1)
% q_method   : method to use for FM estimation (1: CRE1, 2(default):CRE2t, 3: CRE2w 4: CRE2r)
% if_method  : method to use for IF estimation (1: default Biased estimator, 2: unbiased complex estimator)
% q_threshold: threshold applied on window values (default: 10^(-4))
%
%
% OUTPUT:
% tfr:          recursive STFT
% stfr:         recursive vertical synchrosqueezed STFT
% nfreqs:       vector of normalized frequencies (mi:mf)/M
% below, over:  stored out of bounds reassigned energy
% q_hat_map:    continuous time estimated CRE (\hat q)
% if_map:       continuous time estimated IF (\hat\omega)
%
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 27-06-2016
% Ref: [D. Fourer, F. Auger, K.Czarnecki and S. Meignen, Chirp rate and instantaneous frequency estimation. Proc. IEEE ICASSP 2017]


if ~exist('q_method', 'var')
  q_method = 2;
end
if ~exist('if_method', 'var')
 if_method = 1;
end
if ~exist('q_threshold', 'var')
 %% CRE threshold
 q_threshold = 1e-4;  %eps
end

order3 = false;

N        = length(x);
Nfbins   = mf-mi+1;
nfreqs   = (mi:mf)/M;   %% real normalized frequency
ma       = mi:mf;       %% absolute frequency index

q_hat_map= zeros(Nfbins, N);
if_map   = zeros(Nfbins, N);
tfr      = zeros(Nfbins, N);
stfr     = zeros(Nfbins, N);
tfr_k    = zeros(Nfbins, N);
tfr_kp1  = zeros(Nfbins, N);
tfr_kp2  = zeros(Nfbins, N);

if q_method == 4
 order3 = true;
end

if order3
 tfr_kp3  = zeros(Nfbins, N);  
end

if k >=2, tfr_km1  = zeros(Nfbins, N); end
if k >=3, tfr_km2  = zeros(Nfbins, N); end
below    = zeros(1,N);
over     = zeros(1,N);

for m = 1:Nfbins,
 lambda = (mi+m-1)/M;
 pTs    = -1.0/L + 1i*2*pi*lambda;
 alpha  = exp(pTs);
 [a,b]  = Gk2(k,   L, alpha); tfr_k(m,:)   = filter(b,a,x);
 [a,b]  = Gk2(k+1, L, alpha); tfr_kp1(m,:) = filter(b,a,x);
 [a,b]  = Gk2(k+2, L, alpha); tfr_kp2(m,:) = filter(b,a,x);
 if order3
  [a,b]  = Gk2(k+3, L, alpha); tfr_kp3(m,:) = filter(b,a,x);
 end
 if k >=2,  [a,b]  = Gk2(k-1, L, alpha); tfr_km1(m,:) = filter(b,a,x); end
 if k >=3,  [a,b]  = Gk2(k-2, L, alpha); tfr_km2(m,:) = filter(b,a,x); end
 
 tfr(m,:) =  tfr_k(m,:);
 for n = 1:N
  if (abs(tfr_k(m,n)) > eps)
   Tsm1_yTg  = k*L*tfr_kp1(m,n);
   Tsm2_yT2g = k*(k+1)*L^2*tfr_kp2(m,n);
   if k >= 2
    Ts_yDg    = tfr_km1(m,n)/L + pTs*tfr_k(m,n);
   elseif k == 1
    Ts_yDg    = pTs*tfr_k(m,n);   
   end
   if k >= 3
    Ts2_yD2g  = tfr_km2(m,n)/L^2 + 2*pTs*tfr_km1(m,n)/L + pTs^2*tfr_k(m,n);
   elseif k == 2
    Ts2_yD2g  = 2*pTs/L * tfr_km1(m,n) + pTs^2*tfr_k(m,n);
   elseif k == 1
     Ts2_yD2g = pTs^2*tfr_k(m,n);  
   end
   yDTg = k*tfr_k(m,n) + k*L*pTs*tfr_kp1(m,n);
   yTDg = yDTg - tfr_k(m,n);    %yTDg = (k-1)*tfr_k(m,n) + k*L*pTs*tfr_kp1(m,n);
   
   if order3
    Tsm3_yT3g   = L^3 * k * (k+1) * (k+2) * tfr_kp3(m,n);    %L * (k+2)*Tsm2_yT2g;
    %Tsm3_yT3g   = L * (k+2)*Tsm2_yT2g;
    Tsm1_yT2Dg  = L * k * (k-1) * tfr_kp1(m,n) + pTs * L^2 * k * (k+1) * tfr_kp2(m,n);
    Tsm1_yDT2g  = 2 * Tsm1_yTg + Tsm1_yT2Dg;
   end
   
   %% used by LM reassignment algorithm
   %Tsm1_dpsi_domega    =  real(Tsm1_yTg /tfr_k(m,n));
   %Ts_dpsi_dt          =  imag(Ts_yDg   /tfr_k(m,n));
   %d2psi_dtdomega      =  real(yDTg     /tfr_k(m,n) - (Ts_yDg * Tsm1_yTg)/ tfr_k(m,n)^2);
   %Ts2_d2psi_dt2       =  imag(Ts2_yD2g /tfr_k(m,n) - (Ts_yDg            / tfr_k(m,n))^2);
   %Tsm2_d2psi_domega2  = -imag(Tsm2_yT2g/tfr_k(m,n) - (Tsm1_yTg          / tfr_k(m,n))^2);
   
   ntilde = n - Tsm1_yTg/tfr_k(m,n);
   nhat   = round(real(ntilde));
   
   %nhat   = n - round(real(Tsm1_yTg/tfr_k(m,n)));
   
   
   %nhat   = n - real(Tsm1_yTg/tfr_k(m,n));
   mrtilde = M*(Ts_yDg/tfr_k(m,n))/(2.0*pi);
   mrhat  = round(imag(mrtilde)) - mi+1;
   
   %mrhat = round(M*imag(Ts_yDg/tfr_k(m,n))/(2.0*pi)) - mi;
   %mrhat  = M*imag(Ts_yDg/tfr_k(m,n))/(2.0*pi) - mi;
   
   
   
   if q_method == 1   % CRE1
    q_hat_denum = imag(Tsm1_yTg * conj(tfr_k(m,n)));
   elseif q_method == 2 % CRE2_t
    q_hat_denum = (Tsm1_yTg * Ts_yDg - (yDTg-tfr_k(m,n)) * tfr_k(m,n) );
   elseif q_method == 3 %CRE2_w
    q_hat_denum = Tsm1_yTg^2 - Tsm2_yT2g * tfr_k(m,n);
   elseif q_method == 4 % CRE_robust
%     b2 = 1j * (Tsm3_yT3g * tfr_k(m,n) - Tsm2_yT2g * Tsm1_yTg) / tfr_k(m,n)^2;
%     a1 = (Tsm1_yTg * Ts_yDg - yTDg * tfr_k(m,n)) / tfr_k(m,n)^2;
%     
%     b1 = (Tsm1_yDT2g * tfr_k(m,n) - Tsm2_yT2g * Ts_yDg) / tfr_k(m,n)^2;
%     %b1 = (Tsm1_yT2Dg * tfr_k(m,n) - Tsm2_yT2g * Ts_yDg) / tfr_k(m,n)^2;
%     a2 = 1j * (Tsm1_yTg^2 - Tsm2_yT2g * tfr_k(m,n)) / tfr_k(m,n)^2;
%     
%     q_hat_denum = b2 * a1 - b1 * a2;
    
    Ar = [Tsm2_yT2g    -Tsm1_yTg       tfr_k(m,n);...
          Tsm1_yDT2g   -yTDg           Ts_yDg;...
          Tsm3_yT3g    -Tsm2_yT2g      Tsm1_yTg];
    %Tsm1_yDT2g Tsm1_yT2Dg
    q_hat_denum = det(Ar);
   end

   if abs(q_hat_denum) > q_threshold   %% check if signal is FM
     if q_method == 1    %CRE1
       alpha_hat = real(Ts_yDg * conj(tfr_k(m,n))) / q_hat_denum;
       q_hat = 1j * alpha_hat; 
       %q_hat = (mrtilde - 1j*mi)/ntilde;
     elseif q_method == 2 %CRE2_t
       %q_hat =  imag((Ts2_yD2g * tfr_k(m,n) - Ts_yDg^2) / q_hat_denum);  %modified
       q_hat =  (Ts2_yD2g * tfr_k(m,n) - (Ts_yDg)^2) / q_hat_denum;  %modified
       alpha_hat = imag(q_hat);
     elseif q_method == 3 %CRE2_w
       q_hat =  (yTDg * tfr_k(m,n) + tfr_k(m,n).^2 - Tsm1_yTg * Ts_yDg) / q_hat_denum;
       alpha_hat = imag(q_hat);
     elseif q_method == 4 %CRE2_robust
       %y1 = (Ts2_yD2g * tfr_k(m,n) - Ts_yDg^2) / tfr_k(m,n)^2;
       %y2 = 1j * ((yTDg *  tfr_k(m,n)) - Ts_yDg * Tsm1_yTg) / tfr_k(m,n)^2;
       
       u = [Ts_yDg Ts2_yD2g yTDg+tfr_k(m,n)].';
       Ar_inv = pinv(Ar);
       rx = Ar_inv(1,:) * u;
       q_hat = Ar_inv(2,:) * u - 2 * rx * (n-ntilde);
       if abs(q_hat) > 1/q_threshold
         q_hat = 0;
       end
       %q_hat = (b2*y1 - b1 * y2) / q_hat_denum;
       alpha_hat = imag(q_hat);
     end
   else
     q_hat     = 0;
     alpha_hat = 0;
   end
   q_hat_map(m,n) = q_hat;

   %% IF 1, biased when Ax is not constant
   if if_method == 1
     %% normalized frequency
     if_map(m,n) = (imag(mrtilde) + M / (2*pi) * alpha_hat * real(n-ntilde)-mi)/M;
     mrhat_q     = round(imag(mrtilde) + M / (2*pi) * alpha_hat * real(n-ntilde))-mi+1;
     %%mrhat_q = round(mrhat + M / (2*pi) * q_hat * (n-nhat)) 
   
   %% IF 2, unbiased
   elseif if_method == 2 
     %mrhat_q = imag(mrtilde + round(M / (2*pi) * q_hat * (n-ntilde)))-mi;
     
     if_map(m,n) = (imag(mrtilde + M / (2*pi) * q_hat * (n-ntilde))-mi)/M;
     mrhat_q     = round(imag(mrtilde + M / (2*pi) * q_hat * (n-ntilde)))-mi+1;
     %mrhat_q = mrhat + round(imag(M / (2*pi) * q_hat) *real(n-ntilde)); %-mi;
   elseif if_method == 3
     %% IF 3, (classical IF (reassignment operator))
     mrhat_q     = mrhat;
     if_map(m,n) = (imag(mrtilde)-mi)/M;
   end
       

   %denom = (mu + d2psi_dtdomega) * (1 + mu - d2psi_dtdomega) + Ts2_d2psi_dt2 * Tsm2_d2psi_domega2;
   %nhat = n - round((Tsm1_dpsi_domega * (1 + mu - d2psi_dtdomega) - (Tsm2_d2psi_domega2) * (2 * pi * lambda- Ts_dpsi_dt))/denom );

   %Ts_omega_tilde = ((mu^2-d2psi_dtdomega^2) + Ts2_d2psi_dt2 * Tsm2_d2psi_domega2)*2*pi*lambda/denom;
   %Ts_omega_tilde = Ts_omega_tilde - (Tsm1_dpsi_domega * Ts2_d2psi_dt2  - (mu + d2psi_dtdomega)*Ts_dpsi_dt)/denom;
   %mrhat =  round(M*Ts_omega_tilde/(2.0*pi))-mi+1; % eq. 25

%    %% correct m bounds
%    if (mi+mrhat_q>M), mrhat_q=mrhat_q-M; end;
%    if (mrhat_q<0), mrhat_q=0; end;
%    if (mrhat_q>Nfbins), mrhat_q=Nfbins+1; end;
      %% correct m bounds
   %if (mrhat_q>M), mrhat_q=mrhat_q-M; end;
   if (mrhat_q<0), mrhat_q=1; end;
   if (mrhat_q>Nfbins), mrhat_q=Nfbins; end;
   
   %% correct n bounds
   %if (nhat<0), nhat=0;   end;
   %if (nhat>N), nhat=N+1; end;
   if (mrhat_q<=0)
     below(n) = below(n) + tfr(m,n);
   elseif (mrhat_q>=Nfbins+1)
     over(n) = over(n) + tfr(m,n);
   else
     %rtfr(mrhat,nhat) = rtfr(mrhat,nhat) + tfr(m,n);
     stfr(mrhat_q, n) = stfr(mrhat_q, n) + tfr(m,n) * exp(-2*pi*1i*ma(m)/M * n0);
   end;
  end; 
 end;
end;

end

