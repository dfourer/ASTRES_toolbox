function [P1,P2] = proj_T_opt_parallel_xy(t1,t2,x,indT1,indT2) %projection on centro-symmetric toeplitz matrix space
N=size(x,1);
numpatchx=size(t1,3);
numpatchy=size(t1,4);


P1tmp=zeros(6*N,N-2);
P2tmp=zeros(6*N,N-2);
P1=zeros(6*N,N-2,numpatchx,numpatchy);
P2=zeros(6*N,N-2,numpatchx,numpatchy);

for k=1:numpatchx
    for l=1:numpatchy
    for i=1:N^2
        indiceT1=indT1(:,i);
        indiceT2=indT2(:,i);
        indiceT1=indiceT1(find(indiceT1>0));
        indiceT2=indiceT2(find(indiceT2>0));
        t1kl=t1(:,:,k,l);
        t2kl=t2(:,:,k,l);
        P1tmp(indiceT1) = mean([t1kl(indiceT1);t2kl(indiceT2)]);
        P2tmp(indiceT2) = mean([t1kl(indiceT1);t2kl(indiceT2)]);
        P1(:,:,k,l)=P1tmp(:,:);
        P2(:,:,k,l)=P2tmp(:,:);
    end
    end
end

end

% t1h = t1(1:3*N,:);
% t1b = t1(3*N+1:end,:);
% t1b = rot90(t1b,2);
% t2h = t2(1:3*N,:);
% t2b = t2(3*N+1:end,:);
% t2b = rot90(t2b,2);
% P1=zeros(3*N,N-2);
% P2=zeros(3*N,N-2);
% %projection of t1
% for i=1:N
%     P1(3+3*(i-1),1) = (t1h(3+3*(i-1),1) + t1b(3+3*(i-1),1) + t2h(
%     
%     
%     
%     
%     P1(2+3*(i-1),1) = (t1(2+3*(i-1),1)+t1(3+3*(i-1),2)+t1());
% end
% for i=1:2*N
%     P1(2+3*(i-1),1)=(t1(2+3*(i-1),1)+t1(3+3*(i-1),2))/2;
%     P1(3+3*(i-1),2)=(t1(2+3*(i-1),1)+t1(3+3*(i-1),2))/2;
%     P1(1+3*(i-1),N-3)=(t1(1+3*(i-1),N-3)+t1(2+3*(i-1),N-2))/2;
%     P1(2+3*(i-1),N-2)=(t1(1+3*(i-1),N-3)+t1(2+3*(i-1),N-2))/2;
%     for j = 1:N-4
%         P1(1+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%         P1(2+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%         P1(3+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%     end
% end
% %projection of t2
% P2=t2;
% for i=1:2*N
%     P2(2+3*(i-1),1)=(t1(2+3*(i-1),1)+t1(3+3*(i-1),2))/2;
%     P2(3+3*(i-1),2)=(t1(2+3*(i-1),1)+t1(3+3*(i-1),2))/2;
%     P2(1+3*(i-1),N-3)=(t1(1+3*(i-1),N-3)+t1(2+3*(i-1),N-2))/2;
%     P2(2+3*(i-1),N-2)=(t1(1+3*(i-1),N-3)+t1(2+3*(i-1),N-2))/2;
%     for j = 1:N-4
%         P2(1+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%         P2(2+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%         P2(3+3*(i-1),j)=(t1(1+3*(i-1),j)+t1(2+3*(i-1),j+1)+t1(3+3*(i-1),j+2))/3;
%     end
%     P2= [P2 ; rot90(P2,2)];
% end
% end