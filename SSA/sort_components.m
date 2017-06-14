function [ Y, I ] = sort_components( Y )
% [ Y ] = sort_components( Y )
%
% sort component by descending order of frequency

if size(Y,1) < size(Y,2)
  Y = Y.';
end
    
[N, nc] = size(Y);
Nh = round(N/2);
max_f = zeros(1, nc);

for i = 1:nc
 tmp = abs(fft(Y(:,i)));
 [~, max_f(i)] = max(tmp(1:Nh));
end

[~, I] = sort(max_f, 'descend');

Y = Y(:,I);

end

