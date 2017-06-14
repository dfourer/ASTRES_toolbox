%function [A,normA, alpha,beta, ind] = defineParameterMax_line(image)
function [A,normA] = defineParameterMax_diagonal1(image)


[Nx,Ny] = size(image);
N=Nx;
A = spalloc(Nx*Ny,Nx*Ny,4*Nx*Ny);
notnull = 0;

%Imax=zeros(Nx,Ny);
%Imin=zeros(Nx,Ny);

for k=-(Nx-1)+2:(Ny-1)-2
    diagonal = transpose(diag(image,k));
    [indmin, indmax, indzer] = extr(diagonal);
    if (numel(indmin) ~= 0 & numel(indmax) ~= 0)
    ind = sort([indmin,indmax]);
    %alpha = zeros(1,length(ind));
    %beta = zeros(1,length(ind));
    d = 2*ones(1,length(ind));
    d1 = ones(1,length(ind));
    d2 = ones(1,length(ind));
    j=2;
    Ndiag = numel(diagonal);
    %if numel(ind) ~= 0
    notnull = 1;
    for i=1:Ndiag
        if (ind(j) == i)    
            d(j)  = ind(j+1) - ind(j-1);
            d1(j) = ind(j)   - ind(j-1);
            d2(j) = ind(j+1) - ind(j);
      
%             alpha(j) = x(ind(j));
%             beta(j) = (d1(j)*line(ind(j+1)) + d2(j)*line(ind(j-1)))/d(j); 
%             if j==2
%                 alpha(j-1) = alpha(j);
%                 beta(j-1)  = beta(j);   
%             elseif j==length(ind)-1
%                 alpha(j+1) = alpha(j);
%                 beta(j+1)  = beta(j);           
%             end
            if (k < 0)
                A(Nx*(i-1)-k+i,Nx*(ind(j)-1)-k+ind(j)) = 1;
                A(Nx*(i-1)-k+i,Nx*(ind(j-1)-1)-k+ind(j-1)) = d1(j)/d(j);
                A(Nx*(i-1)-k+i,Nx*(ind(j+1)-1)-k+ind(j+1)) = d2(j)/d(j);
                %A(-k+i,ind(j)) = 1;
                %A(-k+i,ind(j-1)) = d1(j)/d(j);
                %A(-k+i,ind(j+1)) = d2(j)/d(j);
            else A(Nx*(k+i-1)+i,Nx*(k+ind(j)-1)+ind(j)) = 1;
                A(Nx*(k+i-1)+i,Nx*(k+ind(j-1)-1)+ind(j-1)) = d1(j)/d(j);
                A(Nx*(k+i-1)+i,Nx*(k+ind(j+1)-1)+ind(j+1)) = d2(j)/d(j);
                %A(i,k+ind(j)) = 1;
                %A(i,k+ind(j-1)) = d1(j)/d(j);
                %A(i,k+ind(j+1)) = d2(j)/d(j);
            end
            %A(Nx*(i-1)+k,Nx*(ind(j)-1)+k) = 1;
            %A(Nx*(i-1)+k,Nx*(ind(j-1)-1)+k) = d1(j)/d(j);
            %A(Nx*(i-1)+k,Nx*(ind(j+1)-1)+k) = d2(j)/d(j);
            j=j+1;
        end
        if j == length(ind);
            break;
        end
    end
    end
end

normA=0;
xn = 10*rand(Nx*Ny,1);
rho_n=1 + 1e-3; 
rho_n1=1;
while abs(rho_n1 - rho_n)>1e-4*rho_n
    rho_n = rho_n1;
    xn1 = A' *A * xn;
    rho_n1 = norm(xn1)/norm(xn);
    %fprintf('normA = %3.2f\n',rho_n1);
    xn = xn1;
end
if rho_n1>normA
    normA=rho_n1;
end

%Normalisation of matrix A
A(:)=A(:)/sqrt(normA);

%save('matriceDemd.mat','A','normA','alpha','beta','ind');
save('matriceDdiagonal1.mat','A','normA');