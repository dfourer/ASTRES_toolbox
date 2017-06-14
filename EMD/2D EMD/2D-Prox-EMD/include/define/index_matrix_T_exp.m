function [indT1,indT2] = index_matrix_T_exp(x)
N=size(x,1);
index=zeros(N,N);

for i =1 : N^2
    index(i)=i;
end

indexT1=operatorT1_exp(index);
indexT2=operatorT2_exp(index);
indT1=zeros(3,N^2);
indT2=zeros(3,N^2);
ind1=zeros(3,1);
ind2=zeros(3,1);

for i=1:N^2
    ind1(:)=0;
    ind2(:)=0;
    tmp1=find(indexT1 == i);
    tmp2=find(indexT2 == i);
    ind1(1:size(tmp1))=tmp1(1:size(tmp1));
    ind2(1:size(tmp2))=tmp2(1:size(tmp2));
    indT1(:,i)=ind1(:);
    indT2(:,i)=ind2(:);
    %indT1(:,i) = find(indexT1 == i);
    %indT2(:,i) = find(indexT2 == i);
end


end