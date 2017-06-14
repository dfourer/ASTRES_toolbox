function imf = synth_tx(Tx,dt,mywav,nv,Cs,clwin)
% synth_mult : extracts ridges contained in Cs on SQ transform Tx and
% computes the modes.

[nr N] = size(Cs);
na = size(Tx,1);
imf = zeros(nr,N);

for j=1:nr
    tmp = zeros(size(Tx));
    for b=1:N
        tmp(max(1,Cs(j,b)-clwin):min(na,Cs(j,b)+clwin),b) = Tx(max(1,Cs(j,b)-clwin):min(na,Cs(j,b)+clwin),b);
    end
    
    imf(j,:) = synth_sq(tmp,dt,mywav,nv);
end