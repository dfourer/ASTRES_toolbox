function [] = eps2pdf( chemin, force )
%function [] = eps2pdf( chemin )
%
% convert all .eps in given folder to .pdf
% (Requires ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: 12/10/2015

if ~exist('force', 'var')
 force = 0;
end

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
  dest_folder     = strcat(chemin,'/pdf/');
  if ~exist(dest_folder, 'dir')
   mkdir(dest_folder);
  end
  target_filename = strcat(dest_folder,filename,'.pdf');

  if exist(target_filename, 'file') && ~force               %% already exist
    fprintf(1,  'file %s already exists\n', target_filename);
    continue;
  end
  
  %% conversion
  fprintf(1,  'generating file %s, ', target_filename);
  com = sprintf('%s %s "%s" "%s"', command, param, source_filename, target_filename);
  [err, res] = system(com);
%   if err ~= 1 && isunix || err ~= 0 && ~isunix
%    error(sprintf('%s\n%s\n', com, res)); 
%   end
end
fprintf('\n');
end
