function [y,d]=ssa(x,L,I)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------                           
%    Author: Francisco Javier Alonso Sanchez    e-mail:fjas@unex.es
%    Departament of Electronics and Electromecanical Engineering
%    Industrial Engineering School
%    University of Extremadura
%    Badajoz
%    Spain
% -----------------------------------------------------------------
%
% SSA generates a trayectory matrix X from the original series x1
% by sliding a window of length L. The trayectory matrix is aproximated 
% using Singular Value Decomposition. The last step reconstructs
% the series from the aproximated trayectory matrix. The SSA applications
% include smoothing, filtering, and trend extraction.
% The algorithm used is described in detail in: Golyandina, N., Nekrutkin, 
% V., Zhigljavsky, A., 2001. Analysis of Time Series Structure - SSA and 
% Related Techniques. Chapman & Hall/CR.

% x1 Original time series (column vector form)
% L  Window length
% y  Reconstructed time series
% r  Residual time series r=x1-y
% vr Relative value of the norm of the approximated trajectory matrix with respect
%	 to the original trajectory matrix

% The program output is the Singular Spectrum of x1 (must be a column vector),
% using a window length L. You must choose the components be used to reconstruct 
% the series in the form [i1,i2:ik,...,iL], based on the Singular Spectrum appearance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Build trajectory matrix

   N=length(x); 
   if L>N/2; L=N-L; end
   K=N-L+1; 
   X=zeros(L,K);                 % [L,K]=size(X)
   for i=1:K
	X(1:L,i)=x(i:L+i-1); 
   end;
    
% Step 2: SVD

   S=X*X'; 
   [EigenVectors,EigenValues]=eig(S);
     
   [d,i]=sort(-diag(EigenValues));  % sort(X) sorts the elements of X in ascending order.
   d=-d; EigenVectors=EigenVectors(:,i); 

%    
   %sev=sum(d);   
   %figure(1); plot(1:L, (d/sev)*100, 1:L, (d/sev)*100, 'rx');
   %title('Singular Spectrum');
   %xlabel('Eigenvalue Number'); ylabel('Eigenvalue (% Norm of trajectory matrix retained)');
   
   V=(X')*EigenVectors;
   % rc=EigenVectors*V';

% Step 3: Grouping

   % I=input('Choose the agrupation of components to reconstruct the series in the form I=[i1,i2:ik,...,iL]');
   Vt=V';
   rca=EigenVectors(:,I)*Vt(I,:);

% Step 4: Reconstruction

   y=zeros(N,1);
   Lp=min(L,K);
   Kp=max(L,K);

   for k=1:Lp-1
     for m=1:k
      y(k)=y(k)+rca(m,k-m+1)/k;
     end
   end

   for k=Lp:Kp
     for m=1:Lp
      y(k)=y(k)+rca(m,k-m+1)/Lp;
     end
   end

   for k=Kp+1:N
     for m=k+1-Kp:N+1-Kp
      y(k)=y(k)+rca(m,k+1-m)/(N-k+1);
     end
   end
end
