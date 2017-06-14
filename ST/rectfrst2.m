function [s_hat] = rectfrst2(tfr, w0T, mm, M)
% [s_hat] = rectfrst2(stfr, w0T, T)
%
% Signal reconstruction from the Stockwellogram (method 2)
% uses the simplified reconstruction formula
% 
% INPUT:
% stfr   : input synchrosqueezed Stockwellogram
% w0T    : central frequency * window duration parameter:  w0 * T, (default: 2pi 50)
% mm     : frequency bin to reconstruct in [-M/2, M/2]
% M      : number of frequency bins
%
% OUTPUT:
% s_hat : reconstructed signal
%
% Author: D.Fourer
% Date: 05-06-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. Signal Processing. 2015]

N = size(tfr,2);

if ~exist('mm', 'var')
  M = size(tfr,1);
  mm = m_axis(M);
end

%% remove zeros (both methods seem to obtain equivalent results)
% method 1 (add epsilon to avoid divide by zero)
%mm = mm + eps;

% method 2 (remove zero from indices, more secure)
I_nz = find(abs(mm) > eps);
mm = mm(I_nz);
tfr = tfr(I_nz, :);
mm = transpose(mm);


R = zeros(1,N);
for n = 1:N 
  R(n) = sum(1./abs(mm) .* exp(1i * 2*pi*(n-1) * mm / M) .* tfr(:, n));
end

%% 2 Compute Ch
dK = 0.01;
K  = 100;
k =  ((w0T-K):dK:(w0T+K))+eps;
C_h = dK * sum( 1./abs(k) .* exp(-(k - w0T).^2 /2));

s_hat = 1 / C_h * R ;

end