% tests if there are enough (3) extrema to continue the decomposition
% Written by Gabriel Rilling
function stop = stop_EMD(r)
[indmin,indmax] = extr(r);
ner = length(indmin) + length(indmax);
stop = ner < 3;
end
