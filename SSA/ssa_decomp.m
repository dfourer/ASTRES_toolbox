function modes = ssa_decomp(x, L, nc, epsilon, force_supk)
% function modes = ssa_decomp(x, L, nc, epsilon, force_supk)
%
% SSA with hierarchical clustering  for components grouping
% In the hierarchical clustering, the last components whose energy contribution is relatively
% null are eliminated from the clustering procedure. epsilon serves as threshold. 
%
% Input:
% x: input signal
% L: SSA window length
% nc: number of components to extract
% epsilon: threshold applied on the values in the singular spectrum (default: 3e-2)
% force_supk: constraint the number of singular values to use (default: supk depends on epsilon)
%
%
% Authors: J. Harmouche and D. Fourer
% Date: Apr. 2016
%
% Refs: 
% [J. Harmouche, D. Fourer, P. Flandrin, F. Auger and P. Borgnat. One
% or Two Components ? The Singular Spectrum Analysis answers. Proc. SLRA'2015. Grenoble, France.June 2015]
% [J. Harmouche, D. Fourer, P. Flandrin, F. Auger and P. Borgnat. Une ou deux composantes:
% la réponse de l'analyse spectrale singulière. Proc. GRETSI'15. Lyon, France]

if ~exist('epsilon', 'var')
 epsilon = 3e-3; %0.5; %
end

if ~exist('force_supk', 'var')
 force_supk = inf;    
end


Nx=length(x);

modes = zeros(nc, 1);
if nc < 1
 return;
end

[~,d]=ssa(x,L,1);

% Remove negligible components to improve clustering
supk = length(find( (d/max(d))>epsilon));

if supk <  nc
 warning('epsilon is too high, trying with eps=%.2f', d(nc)/max(d))
 supk = nc;
end

if force_supk < inf
 supk = force_supk;
end

y = zeros(Nx,supk);
for I = 1:supk;
 [y(:,I),~]=ssa(x,L,I); 
end

wcorr = wCorrMat2(x,L,supk);

if ~isempty(wcorr)
 Z = linkage(wcorr);                       
 T = cluster(Z,'maxclust', nc); 
 modes = zeros(Nx, nc);
 for i = 1:nc       
  idx = find(T == i);

  if length(idx) > 1
   modes(:, i) = sum(y(:, idx).');
  else
   modes(:, i) = y(:, idx).';
  end
 end
end

end

function wCorr = wCorrMat2(x,L,supk)

N=length(x);
w=zeros(N,1);
wCorr=zeros(supk,supk);

for i=1 : supk
 [Y(:,i),d] = ssa(x,L,i); 
end  
for i=1:N
 w(i)=min([i L N-i+1]);
end
for i =1:supk
  for j= 1 : supk
    wCorr(i,j)=  sum(w.*Y(:,i).*Y(:,j))/(sqrt(sum(w.*Y(:,i).^2)*sum(w.*Y(:,j).^2)));
  end
end
end



