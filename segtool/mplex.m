function res = mplex(center,Nf,Delta,as)    
% res : builds the Time-scale segmentation.

na = length(as);

conflicts = zeros(Nf-1,1);
% check for neighbor conflicts (if 2 modes are frequency too close)
for k=1:Nf-1
    if center(k)+Delta >= center(k+1)-Delta
        conflicts(k) = 1;
    end
end

res = zeros(na,1);
% Interior segmentation
if Nf==1
    res(max(1,center(1)-Delta):min(na,center(1)+Delta)) = 1;
else
    for k=1:Nf-1
        if conflicts(k)
            delt = floor(0.5*(center(k+1)-center(k)));
            res(center(k):center(k)+delt) = k;
            res((center(k+1)-delt:center(k+1))) = k+1;        
        else
            res(center(k):center(k)+Delta) = k;
            res(center(k+1)-Delta:center(k+1)) = k+1;
        end
    end    
    % Boundaries
  	res(center(Nf):min(na,center(Nf)+Delta)) = Nf;
    res(max(1,center(1)-Delta):center(1)) = 1;
end