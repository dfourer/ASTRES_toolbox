function [ Y ] = sst_comp_ext(stfr, mask, L, trace)
% [ Y ] = sst_comp_ext(stfr, mask)
%

if ~exist('trace', 'var')
  trace = 0;   
end

nb_comp = size(mask, 3);
N       = size(mask, 2);
M       = size(mask, 1);
Y       = zeros(nb_comp, N);

if size(stfr,1) ~= M || size(stfr, 2) ~= N
  error('Invalid mask')   
end


for i = 1:nb_comp
  tmp_tfr = zeros(M, N);
  tmp_tfr(mask(:,:,i)>0) = stfr(mask(:,:,i)>0);
  
  if trace == 1
   figure(100)
   subplot(121)
   imagesc(abs(stfr).^2)
   title('Synchrosqueezed transform')
   subplot(122)
   imagesc(abs(tmp_tfr).^2)
   title(sprintf('Component %d masked TFR', i))
   pause
  end
  
  Y(i, :) = my_rectfrsgab(tmp_tfr, L, M); 
end

end

