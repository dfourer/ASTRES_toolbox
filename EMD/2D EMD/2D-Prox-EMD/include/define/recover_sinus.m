function c=recover_sinus(t)

Nt = size(t,1);
N=Nt/6;
c=zeros(N,N);

index=zeros(N,N);
for i =1 : N^2
    index(i)=i;
end
indexT1=operatorT1(index);

for i = 1:N^2
    indT1 = find(indexT1 == i);
    c(i) = mean(t(indT1));
end

end