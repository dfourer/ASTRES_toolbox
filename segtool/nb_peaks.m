function n = nb_peaks(gamma,Wx,b)
% nbpeaks : number of peaks of the wavelet transform Wx at time b with threshold gamma.
%   implements solution of Meignen, Oberlin and McLaughlin 2012


[na nb] = size(Wx);

if nargin>2
    n = 0;
    for k=1:na-1
        if abs(Wx(k,b))<=gamma && abs(Wx(k+1,b))>gamma
            n = n+1;
        end
    end
else
    n = zeros(nb,1);
    for b=1:nb
        tmp=0;
        for k=1:na-1
            if abs(Wx(k,b))<=gamma && abs(Wx(k+1,b))>gamma
                tmp = tmp+1;
            end
        end
        n(b) = tmp;
    end
end


