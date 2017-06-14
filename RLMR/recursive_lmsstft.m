function [tfr, stfr, nfreqs, below, over] = recursive_lmsstft(x, k, L, mi, mf, M, mu, n0)
% [tfr, stfr, nfreqs, below, over] = recursive_lmsstft(x, k, L, mi, mf, M, mu, n0)
%
% Compute the recursive Levenberg-Marquardt synchrosqueezed STFT
%
% INPUT:
% x:            input signal
% k:            order
% L:            window length
% mi,mf:        initial and final frequency bin in range [0,M]
% M:            number of frequency bins
% mu:           Levenberg-Marquardt damping parameter
% n0:           time-shift used to allow reconstruction with a causal window (should be greater than 0, default: 1)
%
% OUTPUT:
% tfr:          recursive STFT
% stfr:         recursive LM-synchrosqueezed STFT
% nfreqs:       vector of normalized frequencies (mi:mf)/M
% below, over:  stored out of bounds reassigned energy
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

N        = length(x);
Nfbins   = mf-mi+1;
nfreqs   = (mi:mf)/M;
ma       = mi:mf;       %% absolute frequency index

tfr      = zeros(Nfbins, N);
stfr     = zeros(Nfbins, N);
tfr_k    = zeros(Nfbins, N);
tfr_kp1  = zeros(Nfbins, N);
tfr_kp2  = zeros(Nfbins, N);
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
   
   Tsm1_dpsi_domega    =  real(Tsm1_yTg /tfr_k(m,n));
   Ts_dpsi_dt          =  imag(Ts_yDg   /tfr_k(m,n));
   d2psi_dtdomega      =  real(yDTg     /tfr_k(m,n) - (Ts_yDg * Tsm1_yTg)/ tfr_k(m,n)^2);
   Ts2_d2psi_dt2       =  imag(Ts2_yD2g /tfr_k(m,n) - (Ts_yDg            / tfr_k(m,n))^2);
   Tsm2_d2psi_domega2  = -imag(Tsm2_yT2g/tfr_k(m,n) - (Tsm1_yTg          / tfr_k(m,n))^2);

   denom = (mu + d2psi_dtdomega) * (1 + mu - d2psi_dtdomega) + Ts2_d2psi_dt2 * Tsm2_d2psi_domega2;
   
   nhat = n; %n - round((Tsm1_dpsi_domega * (1 + mu - d2psi_dtdomega) - (Tsm2_d2psi_domega2) * (2 * pi * lambda- Ts_dpsi_dt))/denom );
   
   
   if (nhat<0), nhat=0;   end;
   if (nhat>N), nhat=N+1; end;
  
   Ts_omega_tilde = ((mu^2-d2psi_dtdomega^2) + Ts2_d2psi_dt2 * Tsm2_d2psi_domega2)*2*pi*lambda/denom;
   Ts_omega_tilde = Ts_omega_tilde - (Tsm1_dpsi_domega * Ts2_d2psi_dt2  - (mu + d2psi_dtdomega)*Ts_dpsi_dt)/denom;
   mrhat =  round(M*Ts_omega_tilde/(2.0*pi))-mi+1; % eq. 25

   if (mi+mrhat>M), mrhat=mrhat-M; end;
   if (mrhat<0), mrhat=0; end;
   if (mrhat>Nfbins), mrhat=Nfbins+1; end;
        
   if (mrhat==0)
     below(nhat) = below(nhat) + tfr(m,n);
   elseif (mrhat==Nfbins+1)
     over(nhat) = over(nhat) + tfr(m,n);
   else
     %rtfr(mrhat,nhat) = rtfr(mrhat,nhat) + tfr(m,n);
     stfr(mrhat,nhat) = stfr(mrhat,nhat) + tfr(m,n) * exp(-2*pi*1i*ma(m)/M * n0);
   end;
  end; 
 end;
end;

end

