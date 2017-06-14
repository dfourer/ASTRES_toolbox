function imcoherency = compute_coherency(immodel,sizepatch);

[sizex,sizey]=size(immodel);
imsv3=zeros(sizex,sizey);

for indi = floor(sizepatch/2):sizepatch:sizex-floor(sizepatch/2)
    for indj = floor(sizepatch/2):sizepatch:sizey-floor(sizepatch/2)
        patch=immodel(indi-floor(sizepatch/2)+1:indi+floor(sizepatch/2),indj-floor(sizepatch/2)+1:indj+floor(sizepatch/2));
        [U,S,V]=svd(patch);
        d=diag(S);
        imsv3(indi-floor(sizepatch/2)+1:indi+floor(sizepatch/2),indj-floor(sizepatch/2)+1:indj+floor(sizepatch/2))=sum(d(1:2));    
    end
end

imcoherency=zeros(sizex,sizey);
imcoherency=imsv3/max(imsv3(:));
