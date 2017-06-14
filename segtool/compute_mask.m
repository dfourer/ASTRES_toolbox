function SG = compute_mask(Wx,fs,Cs,clwin,dt,gamma)
% function compute_mask : for the method of Brevdo, it computes a
% posteriori the mask SG of the ridges with respect to the wavelet
% transform WX.
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

ns = length(fs);
[nmodes N] = size(Cs);

% Candidate instantaneous frequency
Wxx = unwrap(angle(Wx),[],2);
Wxr = [Wxx(:, end-1:end) Wxx Wxx(:, 1:2)];
w2 = -Wxr(:, 5:end);
w2 = w2 + 8*Wxr(:, 4:end-1);
w2 = w2 - 8*Wxr(:, 2:end-3);
w2 = w2 + Wxr(:, 1:end-4);
w2 = w2 / (12*dt*2*pi);
w2(abs(Wx)<gamma)=NaN;
w2(abs(w2)>N/2) = NaN;
w = w2;

SG = zeros(size(Wx));
for k=1:nmodes
    for b=1:N
        SG(fs(max(1,Cs(k,b)-clwin))<= w(:,b) & w(:,b)<=fs(min(ns,Cs(k,b)+clwin)),b) = k;
    end
end

%figure(); imagesc(SG);