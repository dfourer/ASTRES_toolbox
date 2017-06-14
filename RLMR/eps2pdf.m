function [] = eps2pdf( chemin )
%function [] = eps2pdf( chemin )
%
% convert all .eps in given folder to .pdf
% (Requires ghostscript http://ghostscript.com/download/ ps2pdf)
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

command = 'ps2pdf';

if isunix
 param = '-dEmbedAllFonts=true -dSubsetFonts=true -dPDFSETTINGS=/prepress -dEPSCrop -dAutoRotatePages';
else
 param = '-dEmbedAllFonts#true -dSubsetFonts#true -dPDFSETTINGS#/prepress -dEPSCrop -dAutoRotatePages';
end
[err, res] = system(command);

if err ~= 1 && isunix || err ~= 0 && ~isunix
  error('eps2pdf cannot be found, please install ghostscript') 
end

d = dir(chemin);

for i = 1:length(d)
  if d(i).name(1) == '.' || length(d(i).name) < 4
     continue;      
  end
  
  filename =  d(i).name(1:(end-4));
  ext      = d(i).name((end-3):end);
  if ~strcmpi(ext, '.eps')
    continue;
  end
  
  source_filename = strcat(chemin,'/',filename,'.eps');
  target_filename = strcat(chemin,'/',filename,'.pdf');
  
  if exist(target_filename, 'file') %% already exist
    delete(target_filename);
    %continue;
  end
  
  %% conversion
  com = sprintf('%s %s "%s" "%s"', command, param, source_filename, target_filename);
  [err, res] = system(com);
%   if err ~= 1 && isunix || err ~= 0 && ~isunix
%    error(sprintf('%s\n%s\n', com, res)); 
%   end
end


end

