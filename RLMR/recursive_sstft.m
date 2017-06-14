function [tfr, stfr, nfreqs, below, over] = recursive_sstft(x, k, L, mi, mf, M, n0)
% [tfr, stfr, nfreqs, below, over] = recursive_sstft(x, k, L, mi, mf, M, n0)
%
% Compute the recursive synchrosqueezing STFT
%
% INPUT:
% x:            input signal
% k:            order
% L:            window length
% mi, mf:       initial and final frequency bin
% M:            number of frequency bins
% n0:           time-shift used to allow reconstruction with a causal window (should be greater than 0, default: 1)
%
% OUTPUT:
% tfr:          recursive STFT
% stfr:         recursive synchrosqueezed STFT
% nfreqs:       vector of normalized frequencies (mi:mf)/M
% below, over:  store the out of bounds reassigned energy
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if ~exist('n0', 'var')
 n0 = 1;    
end

N        = length(x);
Nfbins   = mf-mi+1;
nfreqs   = (mi:mf)/M;
ma       = mi:mf;       %% absolute frequency index

tfr      = zeros(Nfbins, N);
stfr     = zeros(Nfbins, N);
tfr_k    = zeros(Nfbins, N);
tfr_kp1  = zeros(Nfbins, N);
if k >=2, tfr_km1  = zeros(Nfbins, N); end

below    = zeros(1,N);
over     = zeros(1,N);


for m = 1:Nfbins,
 lambda = (mi+m-1)/M;
 pTs    = -1.0/L + 1i*2*pi*lambda;
 alpha  = exp(pTs);
 [a,b]  = Gk2(k,   L, alpha); tfr_k(m,:)   = filter(b,a,x);
 [a,b]  = Gk2(k+1, L, alpha); tfr_kp1(m,:) = filter(b,a,x);
 if k >=2,  [a,b]  = Gk2(k-1, L, alpha); tfr_km1(m,:) = filter(b,a,x); end
 tfr(m,:) = tfr_k(m,:);
 
 for n = 1:N,
  if abs(tfr_k(m,n)) ~= 0
%     Tsm1_yTg=k*L*tfr_kp1(m,n);
%     nhat = n - round(real(Tsm1_yTg/tfr_k(m,n))); % eq. 26
%     if (nhat<0), nhat=0;   end;
%     if (nhat>N), nhat=N+1; end;
   nhat = n;        %% no time reassignment for synchrosqueezing
   
   if k >= 2
    Ts_yDg = tfr_km1(m,n)/L + pTs*tfr_k(m,n);
   elseif k == 1
    Ts_yDg = pTs*tfr_k(m,n);   
   end
   mrhat  = round(M*imag(Ts_yDg/tfr_k(m,n))/(2.0*pi)) - mi+1; %% discrete IF fixed
   
   %if (mi+mrhat>M), mrhat=mrhat-M; end;
   if (mrhat<0), mrhat=1; end;
   if (mrhat>Nfbins), mrhat=Nfbins; end; %+1
   if (nhat==0) || (nhat==N+1)
     continue;  
   end
   
   if (mrhat==0)
     below(n)      = below(n)      + abs(tfr(m,n))^2;
   elseif (mrhat==Nfbins+1)
     over(n)       = over(n)       +  abs(tfr(m,n))^2;
   else
     stfr(mrhat,nhat) = stfr(mrhat,nhat) + tfr(m,n) * exp(-2*pi*1i*ma(m)/M * n0);
   end
  end
 end
end

end

