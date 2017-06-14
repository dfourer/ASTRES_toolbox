function SG = compute_maskWx(Wx,Cs,Delta,cpsi,fs)
% function compute_maskWx : computes the mask for the new reconstruction
% not based on w(a,b).
% If hi is the i-th ridge, hi is the inverse transform of Wx * 1(SG==i).
% Only needed for figure scripts, not for the main soft segtool.
%
% Inputs:
%   Wx : Wavelet transform
%   fs : frequency vector
%   Cs : indexes of ridges
%   clwin : clearing window (see Brevdo et al)
%   dt : sample period
%   gamma : threshold for synchrosqueezing
% Outputs:
%   SG : the masking matrice, same size as Wx. The coefficients such that
%   SG==k are the coefficients of the k-th ridge.
% 

[na N] = size(Wx);
nmodes = size(Cs,1);

SG = zeros(na,N);

for k=1:nmodes
    %a = floor(0.5+log(cpsi ./ fs(Cs(k,:)))/log(2));
    for b=1:N
        a = floor(Cs(k,b));%-3*Delta;
        SG(max(1,a-Delta):min(na,a+Delta),b) = k;
	end
 end


