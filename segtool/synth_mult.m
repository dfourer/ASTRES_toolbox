function imf = synth_mult(Wx,dt,mywav,nv,SG,reg,as)
% synth_mult : extracts ridges contained in Cs on SQ transform Tx and
% computes the modes.

[na N] = size(Wx);
nr = max(SG(:));
imf = zeros(nr,N);

for j=1:nr
    tmp = zeros(size(Wx));
    idx = find(abs(SG-j)<=eps);
    tmp(idx) = Wx(idx);
    
    if reg
        % Projection on the RKHS
        W = projrkhs(tmp,mywav,dt,as,N,nv);
    else
        W = tmp;
    end
    imf(j,:) = synth_sq(W,dt,mywav,nv);
end



