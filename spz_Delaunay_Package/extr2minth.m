% Fixed bugs. Now returns indices in the right order
% Rmk: C (and c) designate rows, R (and r) Columns

function [x,y] = extr2minth(M,th)

[C,R] = size(M);

Mid_Mid = zeros(size(M)); % Boolean matrix. True for matrix which min
                          % is at the middle, and max higher than th

for c = 2:C-1   
    for r = 2:R-1
       T = M(c-1:c+1,r-1:r+1) ;
       Mid_Mid(c, r) = ( (min(min(T)) == T(2, 2)) .* (max(max(T))>th));
    end
end

[x, y] = find(Mid_Mid);

