function sigma = estim_sigma(s)
% estim_sigma : estimates the noise deviation, implements solution of
% Meignen, Oberlin and McLaughlin 2012
% Warning, it uses WaveLab850.

try
    qmf = MakeONFilter('Symmlet',4);
    [y sigma] = NormNoise(s(:)',qmf);
catch
    errordlg('Wavelab is not installed on this computer, estimation of noise aborted');
end

