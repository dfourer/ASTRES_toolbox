function [ w ] = gauss_w( t, sigma )
% function [ w ] = gauss_w( t, sigma )
%
% compute a gaussian function

w = 1/(sigma * sqrt(2*pi)) .* exp(-1/2 * (t/sigma).^2);

end
