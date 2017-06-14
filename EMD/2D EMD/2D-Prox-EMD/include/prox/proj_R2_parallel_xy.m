function P = proj_R2_parallel_xy(t) %projection on rank 2 matrix space
[sizex,sizey,numpatchx,numpatchy]=size(t);
P=zeros(sizex,sizey,numpatchx,numpatchy);
for k=1:numpatchx
    for l=1:numpatchy
    [U,S,V]=svd(t(:,:,k,l));
    [NS1,NS2]=size(S);
    for i=3:min(NS1,NS2)
        S(i,i)=0;
    end
    P(:,:,k,l)=U*S*V';
    end
end
end