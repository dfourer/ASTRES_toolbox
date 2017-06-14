function  [A,normA, alpha,beta, ind] = defineParameterMax_2D(image)


[Nx,Ny] = size(image);
N=Nx;
A = spalloc(Nx*Ny,Nx*Ny,4*Nx*Ny);

%if strcmp(ext,'d2');
%    ind = extr2(x);
%else
%    [indmin, indmax, indzer] = extr(x);
%    ind = sort([indmin,indmax]);
%end
%Computing maxima and minima of x

Imax = imregionalmax(image);
Imin = imregionalmin(image);

indmax = find(Imax==1);
%indmax_x = mod(indmax,N);
%indmax_y = 1+floor(indmax./N);

indmin = find(Imin==1);
%indmin_x = mod(indmin,N);
%indmin_y = 1+floor(indmin./N);

ind = sort([indmin;indmax]);



alpha = zeros(1,length(ind));
beta = zeros(1,length(ind));
%d = 2*ones(1,length(ind));
%d1 = ones(1,length(ind));
%d2 = ones(1,length(ind));
%j=2;
param_alpha_max_tab=[];
param_beta_max_tab=[];
param_gamma_max_tab=[];
param_alpha_min_tab=[];
param_beta_min_tab=[];
param_gamma_min_tab=[];

for i=1:size(indmax)
   %i
   indi=indmax(i);
   if mod(indi,N) == 0
       xi=N;
       yi=floor(indi./N);
   else xi=mod(indi,N);
        yi=1+floor(indi./N);
   end
   distance=[];
   for j=1:size(indmin)
       %j
       indj=indmin(j);
       if mod(indj,N) == 0
           xj=N;
           yj=floor(indj./N);
       else xj=mod(indj,N);
            yj=1+floor(indj./N);
       end
       r=sqrt((xi-xj)^2+(yi-yj)^2);
       distance = [distance ; r];
   end
   %indexes of the 3 closest minima
   index1=find(distance == min(distance));
   if numel(index1) > 1
      index1=index1(1); 
   end
   distance(index1)=100000;
   index1=indmin(index1);
   index2=find(distance == min(distance));
   if numel(index2) > 1
      index2=index2(1); 
   end
   distance(index2)=100000;
   index2=indmin(index2);
   %index3=find(distance == min(distance));
   ind3=find(distance == min(distance));
   if numel(ind3) > 1
      ind3=ind3(1); 
   end
   index3=indmin(ind3);
   %weights of the barycentre of the 3 minima
   if mod(index1,N) == 0
       x1=N;
       y1=floor(index1./N);
   else x1=mod(index1,N);
        y1=1+floor(index1./N);
   end
   
   if mod(index2,N) == 0
       x2=N;
       y2=floor(index2./N);
   else x2=mod(index2,N);
        y2=1+floor(index2./N);
   end
   if mod(index3,N) == 0
       x3=N;
       y3=floor(index3./N);
   else x3=mod(index3,N);
        y3=1+floor(index3./N);
   end
   
   %while (x2-x1)/(y2-y1) == (x3-x1)/(y3-y1)
   if (y1~=y2 && y1~=y3 && x1~=x2 && x1~=x3)
   while (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) == 0
       distance(ind3)=10000;
       ind3=find(distance == min(distance));
       if numel(ind3) > 1
           ind3=ind3(1); 
       end
       index3=indmin(ind3);
       if mod(index3,N) == 0
           x3=N;
           y3=floor(index3./N);
       else x3=mod(index3,N);
           y3=1+floor(index3./N);
       end
   end
   end
   [param_alpha,param_beta,param_gamma] = barycentre (xi,yi,x1,y1,x2,y2,x3,y3);
%    param_gamma = 1;
%    if x1 == xi
%         param_beta = param_gamma * ((y3-yi)*(x1-xi) - (x3-xi)*(y1-yi))/((x2-xi)*(y1-yi)-(x1-xi)*(y2-yi));
%         param_alpha = (-param_beta*(y2-yi) + param_gamma*(y3-yi))/(y1-yi);
%    else
%        param_beta = param_gamma * ((x3-xi)*(y1-yi) - (y3-yi)*(x1-xi))/((y2-yi)*(x1-xi)-(y1-yi)*(x2-xi));
%        param_alpha = (-param_beta*(x2-xi) + param_gamma*(x3-xi))/(x1-xi);
%    end
   %param_alpha
   %param_beta
   %param_gamma
   %param_alpha+param_beta+param_gamma
   %Matrix A
   param_alpha_max_tab=[param_alpha_max_tab param_alpha];
   param_beta_max_tab=[param_beta_max_tab param_beta];
   param_gamma_max_tab=[param_gamma_max_tab param_gamma];

   
   A(indi,indi) = 1;
   A(indi,index1) = param_alpha/(param_alpha+param_beta+param_gamma);
   A(indi,index2) = param_beta/(param_alpha+param_beta+param_gamma);
   A(indi,index3) = param_gamma/(param_alpha+param_beta+param_gamma);
   alpha(indi)=image(xi,yi);
   beta(indi)=(param_alpha*image(x1,y1)+param_beta*image(x2,y2)+param_gamma*image(x3,y3))/(param_alpha+param_beta+param_gamma);
end

for i=1:size(indmin)
   %i
   indi=indmin(i);
   if mod(indi,N) == 0
       xi=N;
       yi=floor(indi./N);
   else xi=mod(indi,N);
        yi=1+floor(indi./N);
   end
   distance=[];
   for j=1:size(indmax)
       %j
       indj=indmax(j);
       if mod(indj,N) == 0
           xj=N;
           yj=floor(indj./N);
       else xj=mod(indj,N);
            yj=1+floor(indj./N);
       end
       r=sqrt((xi-xj)^2+(yi-yj)^2);
       distance = [distance ; r];
   end
   %indexes of the 3 closest maxima
   index1=find(distance == min(distance));
   if numel(index1) > 1
      index1=index1(1); 
   end
   distance(index1)=100000;
   index1=indmax(index1);
   index2=find(distance == min(distance));
   if numel(index2) > 1
      index2=index2(1); 
   end
   distance(index2)=100000;
   index2=indmax(index2);
   ind3=find(distance == min(distance));
   if numel(ind3) > 1
      ind3=ind3(1); 
   end
   index3=indmax(ind3);
   %weights of the barycentre of the 3 maxima
   if mod(index1,N) == 0
       x1=N;
       y1=floor(index1./N);
   else x1=mod(index1,N);
        y1=1+floor(index1./N);
   end
   
   if mod(index2,N) == 0
       x2=N;
       y2=floor(index2./N);
   else x2=mod(index2,N);
        y2=1+floor(index2./N);
   end
   if mod(index3,N) == 0
       x3=N;
       y3=floor(index3./N);
   else x3=mod(index3,N);
        y3=1+floor(index3./N);
   end
   if (y1~=y2 && y1~=y3 && x1~=x2 && x1~=x3)
   while (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) == 0
       distance(ind3)=10000;
       ind3=find(distance == min(distance));
       if numel(ind3) > 1
           ind3=ind3(1); 
       end
       index3=indmax(ind3);
       if mod(index3,N) == 0
           x3=N;
           y3=floor(index3./N);
       else x3=mod(index3,N);
           y3=1+floor(index3./N);
       end
   end
   end
   [param_alpha,param_beta,param_gamma] = barycentre (xi,yi,x1,y1,x2,y2,x3,y3);
   param_alpha_min_tab=[param_alpha_min_tab param_alpha];
   param_beta_min_tab=[param_beta_min_tab param_beta];
   param_gamma_min_tab=[param_gamma_min_tab param_gamma];
%    if x1 == xi
%         param_beta = param_gamma * ((y3-yi)*(x1-xi) - (x3-xi)*(y1-yi))/((x2-xi)*(y1-yi)-(x1-xi)*(y2-yi));
%         param_alpha = (-param_beta*(y2-yi) + param_gamma*(y3-yi))/(y1-yi);
%    else
%        param_beta = param_gamma * ((x3-xi)*(y1-yi) - (y3-yi)*(x1-xi))/((y2-yi)*(x1-xi)-(y1-yi)*(x2-xi));
%        param_alpha = (-param_beta*(x2-xi) + param_gamma*(x3-xi))/(x1-xi);
%    end
   %Matrix A
   A(indi,indi) = 1;
   A(indi,index1) = param_alpha/(param_alpha+param_beta+param_gamma);
   A(indi,index2) = param_beta/(param_alpha+param_beta+param_gamma);
   A(indi,index3) = param_gamma/(param_alpha+param_beta+param_gamma);
   alpha(indi)=image(xi,yi);
   beta(indi)=(param_alpha*image(x1,y1)+param_beta*image(x2,y2)+param_gamma*image(x3,y3))/(param_alpha+param_beta+param_gamma);
end

% for i=1:N    
%     if (ind(j) == i)    
%         d(j)  = ind(j+1) - ind(j-1);
%         d1(j) = ind(j)   - ind(j-1);
%         d2(j) = ind(j+1) - ind(j);
%       
%         alpha(j) = x(ind(j));
%         beta(j) = (d1(j)*x(ind(j+1)) + d2(j)*x(ind(j-1)))/d(j); 
%         if j==2
%             alpha(j-1) = alpha(j);
%             beta(j-1)  = beta(j);   
%         elseif j==length(ind)-1
%             alpha(j+1) = alpha(j);
%             beta(j+1)  = beta(j);           
%         end
%         A(i,ind(j)) = 1;
%         A(i,ind(j-1)) = d1(j)/d(j);
%         A(i,ind(j+1)) = d2(j)/d(j);
%         j=j+1;
%     end
%     if j == length(ind);
%         break;
%     end
% end


tmp = sort([alpha;beta]);
alpha = tmp(1,:);
beta = tmp(2,:);
        
        
    

normA=0;
xn = 10*rand(Nx*Ny,1);
rho_n=1 + 1e-3; 
rho_n1=1;
while abs(rho_n1 - rho_n)>1e-4*rho_n
    rho_n = rho_n1;
    xn1 = A' *A * xn;
    rho_n1 = norm(xn1)/norm(xn);
    fprintf('normA = %3.2f\n',rho_n1);
    xn = xn1;
end
if rho_n1>normA
    normA=rho_n1;
end

%Normalisation of matrix A
A(:)=A(:)/sqrt(normA);

%save('matriceDemd.mat','A','normA','alpha','beta','ind');
save('matriceDemd.mat','A','normA','alpha','beta','ind','param_alpha_max_tab','param_beta_max_tab','param_gamma_max_tab','param_alpha_min_tab','param_beta_min_tab','param_gamma_min_tab');

function [param_alpha,param_beta,param_gamma] = barycentre (xi,yi,x1,y1,x2,y2,x3,y3)

% Mx = [xi-x1 xi-x2 xi-x3 ; yi-y1 yi-y2 yi-y3];
% My = [0 ; 0];
% param_vect = linsolve(Mx,My);
% param_alpha = param_vect(1);
% param_beta = param_vect(2);
% param_gamma = param_vect(3);
% param_alpha
% param_beta
% param_gamma

if x1 == x2 && x1 == x3
    if yi == y1
        param_alpha = 1;
        param_beta = 0;
        param_gamma = 0;
    elseif yi == y2
        param_alpha = 0;
        param_beta = 1;
        param_gamma = 0;
    elseif yi == y3
        param_alpha = 0;
        param_beta = 0;
        param_gamma = 1;
    else param_beta = 1;
        param_gamma = 0;
        param_alpha = -(yi-y2)/(yi-y1);
    end
%     if yi ~= y1
%         param_beta = 1;
%         param_gamma = 0;
%         param_alpha = -(yi-y2)/(yi-y1);
%     else  param_gamma = 0;
%          param_alpha = 1;
%          param_beta = -(yi-y1)/(yi-y2);
%     end
elseif y1 == y2 && y1 == y3
    if xi == x1
        param_alpha = 1;
        param_beta = 0;
        param_gamma = 0;
    elseif xi == x2
        param_alpha = 0;
        param_beta = 1;
        param_gamma = 0;
    elseif xi == x3
        param_alpha = 0;
        param_beta = 0;
        param_gamma = 1;
    else param_beta = 1;
        param_gamma = 0;
        param_alpha = -(xi-x2)/(xi-x1);
    end
%     if xi ~= x1
%         param_beta = 1;
%         param_gamma = 0;
%         param_alpha = -(xi-x2)/(xi-x1);
%     else  param_gamma = 0;
%          param_alpha = 1;
%          param_beta = -(xi-x1)/(xi-x2);
%     end
%     elseif yi ~= y2
%         param_gamma = 0;
%         param_alpha = 1;
%         param_beta = -(yi-y1)/(yi-y2);
%     else param_alpha = 1;
%         param_beta = 0;
%         param_gamma = - (yi-y1)/(yi-y3);
%    end
% elseif y1 == y2 && y1 == y3
%     param_beta = 1;
%     param_gamma = 1;
%     if xi ~= x1
%         param_alpha = (x2+x3-2*xi)/(xi-x1);
%     elseif xi ~= x2
%         param_alpha = (x1+x3-2*xi)/(xi-x2);
%     else param_alpha = (x1+x2-2*xi)/(xi-x3);
%     end
elseif (x2-x1)/(y2-y1) == (x3-x1)/(y3-y1) %1,2,3 alignés
    %Eq y=ax+b
    a = (y2-y1)/(x2-x1);
    b = (y1*x2-y2*x1)/(x2-x1);
    %Coordonnées du projeté orthogonal de I:
    xs = (xi*(x2-x1)+(b-yi)*(y2-y1))/(x2-x1+a*(y2-y1));
    ys = a*xs+b;
    if xs == x1
        param_alpha = 1;
        param_beta = 0;
        param_gamma = 0;
    elseif xs == x2
        param_alpha = 0;
        param_beta = 1;
        param_gamma = 0;
    elseif xs == x3
        param_alpha = 0;
        param_beta = 0;
        param_gamma = 1;
    else param_beta = 1;
        param_gamma = 0;
        param_alpha = -(xs-x2)/(xs-x1);
    end          
elseif xi == x1 && xi == x2
    param_gamma = 0;
    param_beta = 1;
    param_alpha = (y2-yi)/(yi-y1);
elseif yi == y1 && yi == y2
    param_gamma = 0;
    param_beta = 1;
    param_alpha = (x2-xi)/(xi-x1);
elseif xi == x1 && xi == x3
    param_beta = 0;
    param_gamma = 1;
    param_alpha = (y3-yi)/(yi-y1);
elseif yi == y1 && yi == y3
    param_beta = 0;
    param_gamma = 1;
    param_alpha = (x3-xi)/(xi-x1);
elseif xi == x2 && xi == x3
    param_alpha = 0;
    param_gamma = 1;
    param_beta = (y3-yi)/(yi-y2);
elseif yi == y2 && yi == y3
    param_alpha = 0;
    param_gamma = 1;
    param_beta = (x3-xi)/(xi-x2);
elseif xi == x1
    param_gamma = 1;
    param_beta = param_gamma * (x3-xi)/(xi-x2);
    param_alpha = (param_beta*(y2-yi) + param_gamma*(y3-yi))/(yi-y1);
elseif yi == y1
    param_gamma = 1;
    param_beta = param_gamma * (y3-yi)/(yi-y2);
    param_alpha = (param_beta*(x2-xi) + param_gamma*(x3-xi))/(xi-x1);
elseif xi == x2
    param_gamma = 1;
    param_alpha = param_gamma * (x3-xi)/(xi-x1);
    param_beta = (param_alpha*(y1-yi) + param_gamma*(y3-yi))/(yi-y2);
elseif yi == y2
    param_gamma = 1;
    param_alpha = param_gamma * (y3-yi)/(yi-y1);
    param_beta = (param_alpha*(x1-xi) + param_gamma*(x3-xi))/(xi-x2);
elseif xi == x3
    param_beta = 1;
    param_alpha = param_beta * (x2-xi)/(xi-x1);
    param_gamma = (param_alpha*(y1-yi) + param_beta*(y2-yi))/(yi-y3);
elseif yi == y3
    param_beta = 1;
    param_alpha = param_beta * (y2-yi)/(yi-y1);
    param_gamma = (param_alpha*(x1-xi) + param_beta*(x2-xi))/(xi-x3);
% elseif (yi-y2)*(xi-x1)+(yi-y1)*(xi-x2) == 0;
%     param_alpha = 1;
%     param_beta = 2;
%     param_gamma = 1;
    %param_gamma = 0;
    %param_beta = 1;
    %param_alpha = (param_beta*(xi-x2) + param_gamma*(xi-x3))/(x1-xi);
     %param_alpha = 1;
     %param_gamma = ((y1-yi)*(xi-x2) + (x1-xi)*(y2-yi))/((x3-xi)*(yi-y2)+(xi-x2)*(yi-y3));
     %param_beta = (x1-xi + param_gamma*(x3-xi))/(xi-x2);
elseif (yi-y2)/(xi-x2) == (yi-y1)/(xi-x1) && (yi-y3)/(xi-x3) == (yi-y1)/(xi-x1)
     param_gamma = 1;
     param_beta = 1;
     param_alpha = - (2*xi-x2-x3)/(xi-x1);
    %param_gamma = 0;
    %param_beta = 1;
    %param_alpha =  - (xi-x2)/(xi-x1);
elseif (yi-y2)/(xi-x2) == (yi-y1)/(xi-x1)
    param_gamma = 0;
    param_beta = 1;
    param_alpha = - (xi-x2)/(xi-x1);
elseif (yi-y3)/(xi-x3) == (yi-y1)/(xi-x1)
    param_gamma = 1;
    param_beta = 0;
    param_alpha = - (xi-x3)/(xi-x1);
elseif (yi-y3)/(xi-x3) == (yi-y2)/(xi-x2)
    param_alpha = 0;
    param_beta = 1;
    param_gamma = - (xi-x2)/(xi-x3);
else param_gamma = 1;
     param_beta = param_gamma * ((xi-x3)*(yi-y1) - (yi-y3)*(xi-x1))/((yi-y2)*(xi-x1)-(yi-y1)*(xi-x2));
     param_alpha = (param_beta*(xi-x2) + param_gamma*(xi-x3))/(x1-xi);
end
