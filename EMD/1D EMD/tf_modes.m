%INPUTS :
%       modes : K x N vector containing the modes
%       M : frequency span
%       method : spectral analysis method
%           'Hilbert' : analytic signal (faster)
%           'Prony' : local Prony spectral analysis (more robust)
%       patchsize : size of patches (only for Prony method)

function [tfr,amp,freq] = tf_modes(modes,M,method,patchsize)

[K,N] = size(modes);

switch nargin
    case 1
        M = N;
        method = 'Hilbert';
        patchsize=10;
    case 2
        method = 'Hilbert';
        patchsize=10;
    case 3
        patchsize=10;
end
amp=zeros(K,N);
freq=zeros(K,N);
frequency=zeros(K,N);
tfr=zeros(M,N);

    
    
if strcmp(method, 'Hilbert')
    modes_a = transpose(hilbert(modes'));
    amp=(abs(modes_a));
    phase=angle(modes_a);
    phase=unwrap(phase);

    for k=1:K
        freq(k,1:N-1)=(diff(phase(k,:))./diff(1:N))/(2*pi);   
    end
end

if strcmp(method,'Prony')
    for k=1:K
        modes_a = transpose(hilbert(modes'));
        amp=(abs(modes_a));
        for ind=1:N-patchsize
            patch=modes(k,ind:ind+patchsize);
            [freq_prony,crit,dc1,dc2,c]=prony_1D(patch,2);
            freq(k,ind+floor(patchsize/2))=abs(freq_prony(1));
        end 
    end
end
    

%compute frequency axe for time-frequency representation
for ind=1:numel(freq)
    if (abs(freq(ind)) >0.5)
        frequency(ind)=1;
    elseif (freq(ind) >= 0)
       frequency(ind) = 1+round(freq(ind)*M);
   else frequency(ind) = M+round(freq(ind)*M);
   end
end


for k=1:K
    for n=1:N-1
        tfr(frequency(k,n),n)=tfr(frequency(k,n),n)+amp(k,n)^2;
    end
end

end

function [freq,crit,dc1,dc2,c]=prony_1D(y,nsin)

N=length(y);

%parameters
mu = 0.1;
gamma = 0.5;
Lmax = 500;

Tdata=operatorT1_1D(y);

%Construction of weight matrix
p1=zeros(3,N-2);
p1(:)=1/6;
p1(1,N-3)=1/4;
p1(1,N-2)=1/2;
p1(2,1)=1/4;
p1(2,N-2)=1/4;
p1(3,1)=1/2;
p1(3,2)=1/4;
p=repmat(p1,2,1);

%criterion
crit = zeros(1,Lmax);
dc1 = zeros(1,Lmax);
dc2 = zeros(1,Lmax);

%initialization
t=zeros(6,N-2);
PT=zeros(6,N-2);
s=zeros(6,N-2);

%Iterations
for l = 1:Lmax
    %l
    
    

    t = proj_rank_1D(s + gamma*(t-s) - mu*p.*(t-Tdata),nsin);
    
    PT = proj_T_1D(2*t-s,y);
    
    s = s - t + PT;
    
    %figure(1);
    %subplot(311);imagesc(t);axis image off;
    %subplot(312);imagesc(PT);axis image off;
    %subplot(313);imagesc(s);axis image off;
    %pause;
    
    %criterion
    M = sqrt(p).*(t-Tdata);
    crit(l)=norm(M,'fro');
    tT=proj_T_1D(t,y);
    dc1(l)=norm(t-tT,'fro');
    dc2(l)=norm(t-proj_rank_1D(t,nsin),'fro')^2;
    if (l>1 & (abs(crit(l)-crit(l-1))<1e-12) & (abs(dc1(l)-dc1(l-1))<1e-12) & (abs(dc2(l)-dc2(l-1))<1e-12))
        break;
    end

end

%Compute sinusoïdal estimate
c=recover_sinus_1D(t);

%Reshape t1 and t2 into 2N2(N1-2)X3 matrix
tr1=zeros(N-nsin,nsin+1);

for indi=1:N-nsin
   tr1(indi,:) = c(N-nsin+1-indi:N+1-indi); 
end

tr2=rot90(tr1,2);
tr=[tr1 ; tr2];

%Compute annihilating filters
[U,S,V]=svd(tr);
h=V(:,nsin+1);

rootsh=roots(h);

freq=angle(rootsh)/(2*pi);

end

function t=operatorT1_1D(x)
N=length(x);
%construction of top matrix
t1=zeros(3,N-2);
for i=1:N
    t1(1,:) = x(3:N);
    t1(2,:) = x(2:N-1);
    t1(3,:) = x(1:N-2);    
    
end
%duplication of top matrix
t2=rot90(t1,2);
%centro-symmetric matrix
t=[t1 ; t2];
end

function P = proj_rank_1D(t,rank) %projection on rank 2 matrix space
[U,S,V]=svd(t);
[NS1,NS2]=size(S);
for i=rank+1:min(NS1,NS2)
    S(i,i)=0;
end
P=U*S*V';
end

function P = proj_T_1D(t,x) %projection on centro-symmetric toeplitz matrix space
N=length(x);
index=zeros(N,1);

for i =1 : N
    index(i)=i;
end

indexT=operatorT1_1D(index);
P=zeros(6,N-2);

for i=1:N
    indT = find(indexT == i);
    P(indT) = mean(t(indT));
end

end

function c=recover_sinus_1D(t)

Nt = size(t,2);
N=Nt+2;
c=zeros(N,1);

index=zeros(N,1);
for i =1 : N
    index(i)=i;
end
indexT=operatorT1_1D(index);

for i = 1:N
    indT = find(indexT == i);
    c(i) = mean(t(indT));
end

end