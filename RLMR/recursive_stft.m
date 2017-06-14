function [tfr, nfreqs] = recursive_stft(x, k, L, mi, mf, M)
% [tfr, nfreqs] = recursive_stft(x, k, L, mi, mf, M)
%
% Computes the recursive STFT
%
% INPUT:
% x:        input signal
% k:        order
% L:        window length
% mi, mf:   initial and final frequency bins
% M:        total number of frequency bins
%
% OUTPUT:
% tfr:      recursive STFT
% nfreqs:   vector of normalized frequencies (mi:mf)/M
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

N      = length(x);
Nfbins = mf-mi+1;
nfreqs = (mi:mf)/M;
tfr    = zeros(Nfbins, N);

for m = 1:Nfbins,
 lambda = (mi+m-1)/M;
 pTs    = -1.0/L + 1i*2*pi*lambda;
 alpha  = exp(pTs);
 [a,b]  = Gk2(k, L, alpha);
 tfr(m,:) = filter(b,a,x);
end;

%tfr = abs(tfr).^2;
end
