function [ text] = gen_tab(vtitle, v_indices, vals, id)
% [ text] = gen_tab(vtitle, v_indices, vals)
% Generate result Table in latex format
% 
% INPUT:
% vtitle: parameter name
% v_indices: vector containing parameter reference values
% v_vals: vector containing the results obtained for each parameter reference value
%
% OUTPUT:
% text: latex code of the resulting tabular
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Real-time recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]
if ~exist('id', 'var')
  id = 'a';
end
 if length(v_indices) ~= length(vals)
  error('v_indices and vals must have the same length')
 end
 
 text = '\\begin{tabular}{|';
 for i = 1:(length(v_indices)+2)
  text =  strcat(text, 'l|');
 end
 text =  strcat(text, '} \\hline\n');

 v_opt = max(vals);

 text =  strcat(text, sprintf('\\\\multirow{2}{*}{\\\\!(%s)\\\\!}\\t& %s\\t\\t ', id, vtitle));
 for i = 1:length(v_indices)
  if abs(fix(v_indices(i))-v_indices(i)) > eps
   v = sprintf('%0.02f', v_indices(i));
  else
   v = sprintf('%d', v_indices(i));
  end

  if vals(i) == v_opt
    text =  strcat(text, sprintf('&\\t\\\\textbf{%s}\\t', v));  
  else
   text =  strcat(text, sprintf('&\t%s\t', v));  
  end
 end
 text =  strcat(text, '\\\\ \n'); %\\hline

 text =  strcat(text, '\t\t\t\t&\\text{RQF}(dB)\t');
 for i = 1:length(vals)
  if abs(fix(vals(i))-vals(i)) > eps
   v = sprintf('%0.02f', vals(i));
  else
   v = sprintf('%d', vals(i));
  end
  
  if vals(i) == v_opt
   text =  strcat(text, sprintf('&\\t\\\\textbf{%s}\\t', v));  
  else
   text =  strcat(text, sprintf('&\t%s\t', v));  
  end
 end
 text =  strcat(text, '\\\\\\hline\n');  
 text =  strcat(text, '\\end{tabular}\n');

end

