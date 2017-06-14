function [s_hat] = rectfrsst(stfr, w0T, mm)
% [s_hat] = rectfrsst(stfr, w0T, T)
%
% Signal reconstruction from the synchrosqueezed S-transform
% 
% INPUT:
% stfr   : input synchrosqueezed Stockwellogram
% w0T    : central frequency * window duration parameter:  w0 * T, (default: 2pi 50)
% mm     : frequency bin to reconstruct in [-M/2, M/2]
%
% OUTPUT:
% s_hat : reconstructed signal
%
% Author: D.Fourer
% Date: 05-06-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. Signal Processing. 2015]


if ~exist('w0T', 'var')
 w0T = 2*pi*50;
end

if ~exist('mm', 'var')
 M = size(stfr, 1);
 mm = m_axis(M);
end

%% remove zeros (both methods seem to obtain equivalent results)
% method 1 (add epsilon to avoid divide by zero)
%mm = mm + eps;

% method 2 (remove zero from indices, more secure)
I_nz = find(abs(mm) > eps);
mm = transpose(mm(I_nz));
stfr = stfr(I_nz, :);

R = sum(stfr .* repmat(1./ abs(mm),1,size(stfr,2)),1);


%% 2 Compute Ch
dK = 0.01;
K  = 100;
k  = ((w0T-K):dK:(w0T+K))+eps;
k  = k(abs(k)>eps);
%plot(k,  1./abs(k) .* exp(-(k - w0T).^2 /2))
%pause
C_h = dK * sum( 1./abs(k) .* exp(-(k - w0T).^2 /2));


s_hat = 1 / (C_h) * R ;
end