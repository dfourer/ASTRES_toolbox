function t=operatorT1(x)
N=size(x,1);
xt=transpose(x);
%construction of top matrix
t1=zeros(3*N,N-2);
for i=1:N
    t1(1+3*(i-1),:) = xt(i,3:N);
    t1(2+3*(i-1),:) = xt(i,2:N-1);
    t1(3+3*(i-1),:) = xt(i,1:N-2);
end
%duplication of top matrix
t2=rot90(t1,2);
%centro-symmetric matrix
t=[t1 ; t2];
end