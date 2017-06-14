function Nf = estim_nf(Wx,as,Gammin,Ndis)
% estim_nf : estimates the number of modes, implements solution of
% Meignen, Oberlin and McLaughlin 2012

[na N] = size(Wx);

ttab = zeros((N)*Ndis,1);
cpt = 1;
for b=1:N
    [vm idx] = max(abs(Wx(:,b)));
    if idx>1 
        Gammax = abs(Wx(idx-1,b));
    else
        Gammax = abs(Wx(idx+1,b));
    end
    for gam = linspace(Gammin,Gammax,Ndis)
        tmp = 0;
        for a=1:na-1
            if abs(Wx(a,b))<gam && abs(Wx(a+1,b))>=gam
                tmp = tmp+1;
            end
        end
        ttab(cpt) = tmp;
        cpt = cpt+1;
    end
end

Nf = floor(median(ttab));

            
            
               
            
            
    



