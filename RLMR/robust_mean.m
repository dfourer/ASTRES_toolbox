function [ x_bar ] = robust_mean( x, alpha )
% [ x_bar ] = robust_mean( x )
%
% return the robust mean excluding outliers
% by considering elements in [c-alpha*s  s+alpha*s]
% where c is the mean of x and s is the standard deviation of x
%

if ~exist('alpha', 'var')
  alpha = 2;   
end

bounds = std(x) * [-alpha alpha] + mean(x);
I = find(and(x >= bounds(1), x <= bounds(2)));

if isempty(I)
 error('Invalid alpha value');
end

x_bar = mean(x(I));

end

