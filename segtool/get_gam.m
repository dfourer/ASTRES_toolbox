function [flag gam] = get_gam(vec,Nf,Ndis,gmin)
% get_gam : find the best threshold for the vector vec in order to
% detecting Nf modes.
% flag = 0 if no convenient threshold exists

na = length(vec);
gvec = linspace(gmin,max(vec),Ndis); 
thresh = [];

for g=gvec
    cpt = 0;
    for a=1:na-1
        if vec(a)<g && vec(a+1)>=g
            cpt = cpt+1;
        end
    end
    if cpt==Nf
        thresh(end+1)=g;
    end
end

if isempty(thresh)
    flag = 0;
    gam = 0;
else
    gam = median(thresh);
    flag = 1;
end
