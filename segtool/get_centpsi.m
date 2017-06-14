function cpsi = get_centpsi(mywav)
% get_centpsi : gives the centre frequency of a wavelet.
% mywav : either a complex Morlet, Generalized Morse or Bump wavelet.

if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    beta = str2num(mywav(v1:v2-2));
    gamma = str2num(mywav(v2:end));
    cpsi = sqrt(beta*gamma);
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    Fc= str2num(mywav(v2:end));
    cpsi = Fc;
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    mu = str2num(mywav(v1:v2-2));
    cpsi = mu;
end
