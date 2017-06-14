function [N x n1] = mypad(s)
% mypad : symmetric padding of vector s.

n = length(s);

% Computing the power of 2
N = 2^(1+round(log2(n+eps)));
n1 = floor((N-n)/2);

% Padding
padtype = 'symmetric';
try
    sleft = padarray(s, n1, padtype, 'pre');
    sright = padarray(s, n1, padtype, 'post');
catch
    warning('Image processing Toolbox is not working. Using constant pading');
    sleft = flipud(s(2:n1+1));
    sright = flipud(s(end-n1:end-1));
end
x = [sleft(1:n1); s; sright(end-n1+1:end)];
