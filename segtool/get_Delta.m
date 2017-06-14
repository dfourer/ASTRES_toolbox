function Delta = get_Delta(mywav)
% get_Delta : gives the bandwidth of a wavelet.
% mywav : either a complex Morlet, Generalized Morse or Bump wavelet.

if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    beta = str2num(mywav(v1:v2-2));
    gamma = str2num(mywav(v2:end));
    warndlg('Not implemented!!!!!');
    Delta = 1;
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    Fb = str2num(mywav(v1:v2-2));
    Delta = 1.6/Fb;
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    sigma= str2num(mywav(v2:end));
    Delta = sigma;
end
