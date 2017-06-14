function [gamma Nf flag] = comp_gamma(Wx,Gam)
% comp_gamma : computes the adaptative threshold, implements solution of
% Meignen, Oberlin and McLaughlin 2012

Gam = sort(Gam);

[na nb] = size(Wx);
nG = length(Gam);

tmpmin = zeros(nG,1);
tmpmean = zeros(nG,1);

for k=1:nG
    n = nb_peaks(Gam(k),Wx);
    tmpmin(k) = min(n);
    tmpmean(k) = mean(n);
end

% the mean value : estimation of the number of modes
Nf = floor(mean(tmpmean)+1-eps);

% Find if ghat exists
if max(tmpmin)<Nf % ghat does not exist. Computing gbar
    flag = 1;
    [k idx] = min(abs(tmpmean-Nf));
    gamma = Gam(idx);
else % ghat exists
    flag = 0;
    tmp = Inf;
    for k=nG:(-1):1
        if tmpmin(k)>=Nf && abs(tmpmean(k)-Nf)<tmp
            tmp = abs(tmpmean(k)-Nf);
            gamma = Gam(k);
        end
    end
end

    
    