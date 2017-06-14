function [amplitude,freq,orientation,coherency,model,crit,dc1,dc2]=prony_2D(imf,sizepatch)

[sizex,sizey]=size(imf);


%size of patch
N=sizepatch;

cinit=zeros(N,N);


%parameters
sigma= 10;
tau = 0.99/(2*sigma + 2);
Lmax = 500;

%number of patches
numpatchx=floor(sizex/N);
numpatchy=floor(sizey/N);
numpatch=numpatchx*numpatchy;

T1data=zeros(6*N,N-2,numpatchx,numpatchy);
T2data=zeros(6*N,N-2,numpatchx,numpatchy);

patch=zeros(N,N);
[indT1,indT2] = index_matrix_T(patch);
indx=1;
indy=1;
for indi = floor(N/2):N:sizex-floor(N/2)
    for indj = floor(N/2):N:sizey-floor(N/2)
        patch=imf(indi-floor(N/2)+1:indi+floor(N/2),indj-floor(N/2)+1:indj+floor(N/2));
        %Compute Toeplitz matrixes of patch
        T1patch=operatorT1(patch);
        T2patch=operatorT2(patch);
        T1data(:,:,indx,indy) = T1patch;
        T2data(:,:,indx,indy) = T2patch;
        %Construction of weight matrix
        p1=zeros(3,N-2);
        p1(:)=1/6;
        p1(1,N-3)=1/4;
        p1(1,N-2)=1/2;
        p1(2,1)=1/4;
        p1(2,N-2)=1/4;
        p1(3,1)=1/2;
        p1(3,2)=1/4;
        p1rep=repmat(p1,2*N,1);
        p=zeros(6*N,N-2,numpatchx,numpatchy);
        for k1=1:numpatchx
            for k2=1:numpatchy
            p(:,:,k1,k2)=p1rep(:,:);
            end
        end
        indy=indy+1;
    end
    indx=indx+1;
    indy=1;
end
%criterion
crit = zeros(1,Lmax);
dc1 = zeros(1,Lmax);
dc2 = zeros(1,Lmax);



%initialization
t1=zeros(6*N,N-2,numpatchx,numpatchy);
t2=zeros(6*N,N-2,numpatchx,numpatchy);
s1=zeros(6*N,N-2,numpatchx,numpatchy);
s2=zeros(6*N,N-2,numpatchx,numpatchy);


%Iterations
for l = 1:Lmax
    t1p = proj_R2_parallel_xy(t1 - 2*tau*p.*(t1-T1data) - tau*s1);
    t2p = proj_R2_parallel_xy(t2 - 2*tau*p.*(t2-T2data) - tau*s2);
    
    %compute projection on T
    [PT1,PT2] = proj_T_opt_parallel_xy(s1+sigma*(2*t1p-t1),s2+sigma*(2*t2p-t2),patch,indT1,indT2);
    
    s1 = s1 + sigma*(2*t1p-t1) - PT1;
    s2 = s2 + sigma*(2*t2p-t2) - PT2;
    
    t1 = t1p;
    t2 = t2p;
    
    %criterion
    M1 = sqrt(p).*(t1-T1data);
    M2 = sqrt(p).*(t2-T2data);
    normM1=zeros(numpatchx,numpatchy);
    normM2=zeros(numpatchx,numpatchy);
    normt1T=zeros(numpatchx,numpatchy);
    normt2T=zeros(numpatchx,numpatchy);
    normt1R2=zeros(numpatchx,numpatchy);
    normt2R2=zeros(numpatchx,numpatchy);
    for k=1:numpatchx
        for k2=numpatchy
        normM1(k1,k2)=norm(M1(:,:,k1,k2),'fro')^2;
        normM2(k1,k2)=norm(M2(:,:,k1,k2),'fro')^2;
        end
    end
    crit(l)=sqrt(sum(normM1(:))+sum(normM2(:)));
    [t1T,t2T]=proj_T_opt_parallel_xy(t1,t2,patch,indT1,indT2);
     for k1=1:numpatchx
        for k2=1:numpatchy
        normt1T(k1,k2)=norm(t1(:,:,k1,k2)-t1T(:,:,k1,k2),'fro')^2;
        normt2T(k1,k2)=norm(t2(:,:,k1,k2)-t2T(:,:,k1,k2),'fro')^2;
        end
    end
    dc1(l)=sqrt(sum(normt1T(:))+sum(normt2T(:)));
    t1R2=proj_R2_parallel_xy(t1);
    t2R2=proj_R2_parallel_xy(t2);
    for k1=1:numpatchx
        for k2=1:numpatchy
        normt1R2(k1,k2)=norm(t1(:,:,k1,k2)-t1R2(:,:,k1,k2),'fro')^2;
        normt2R2(k1,k2)=norm(t2(:,:,k1,k2)-t2R2(:,:,k1,k2),'fro')^2;
        end
    end
    dc2(l)=sqrt(sum(normt1R2(:))+sum(normt2R2(:)));
end

%Compute sinusoidal estimate

[amplitude,orientation,freq,model]=sinus_estimation_parallel(t1,N,sizex,sizey);

coherency = compute_coherency(immodel,sizepatch);


end