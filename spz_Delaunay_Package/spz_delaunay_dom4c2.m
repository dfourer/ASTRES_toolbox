function [ xrpond ] = spz_delaunay_dom4c2( tfr, mask )
%
% P. Flandrin & Ph. Depalle
% 2015, June 10th
%
% inputs
%   tfr = STFT (output from spz_delaunay_dom4a)
%   inseg = masks (output from spz_delaunay_dom4b)
%   domselect = selected sub-domain in inseg (from colorbar)
% output
%   xrpond = reconstructed component

% mask


% masked STFT
[AA,BB] = size(tfr) ;
tfrpond = tfr(1:1+AA/2,:).* mask;

% reconstruction 
xrpond = zeros(1,BB) ;
for n = 1:BB
     xrpond(n) = real(sum(tfrpond(1:1+AA/2,n))) ; % from masked STFT
end

end

