function [ mask ] = ridge_detect_brvmask(stfr, m_axis, nr, lambda, clwin, K)

addpath('segtool');

[M, N] = size(stfr);
[Cs, Es] = brevridge_mult(stfr, m_axis, 0, nr, lambda, clwin);
mask = zeros(M, N, nr);

for i = 1:nr
  
  for n = 1:N    
   if Cs(i,n) > 1 && abs(stfr(Cs(i,n),n)) > lambda
     m_r = unique(max(1,Cs(i,n)-K):min(Cs(i,n)+K,M));
     mask(m_r,n,i) = 1;
   end
  end
end

end
