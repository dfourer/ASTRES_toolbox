function resultsprony = run_prony_2D(imf,sizepatch,K)
%
% The function prony_2D performs a local Prony spectral analysis on each IMF.
% The spectral analysis is performed according to the work 
% "2D Prony-Huang Transform: A New Tool for 2-D Spectral Analysis,
% J. Schmitt, N. Pustelnik, P. Borgnat, P. Fladrin, L.Condat,
% IEEE Transactions on Image Processing, 2014".
%
%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% imf                       : cell of size K containing the IMFs obtained with EMD 
% sizepatch                 : size of patches for each IMF
% K                         : number of IMFs
%
%%%%%%%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% resultsprony              : results of prony spectral analysis
% resultsprony.amplitude    : cell of size K containing the amplitude images of
%                             each IMF
% resultsprony.freq         : frequency images of each IMF
% resultsprony.orientation  : orientation images of each IMF
% resultsprony.coherency    : coherency images of each IMF
% results.model             : denoised IMFs
% results.crit              : convergence criterion of each IMF
% results.dc1               : distance to the set of block-Toeplitz matrix
%                             for each IMF
% results.dc2               : distance to the set of rank 2 matrix
%                             for each IMF
%
%%%%%%%%SUBFONCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% prony_2D : compute local prony spectral analysis on a single IMF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:lambda.K
    [amplitude,freq,orientation,coherency,model,crit,dc1,dc2]=prony_2D(resultsemd.imf{k},sizepatch(k));
    resultsprony.amplitude{k}   = amplitude;
    resultsprony.freq{k}        = freq;
    resultsprony.orientation{k} = orientation;
    resultsprony.coherency{k}   = coherency;
    resultsprony.model{k}       = model;
    resultsprony.crit{k}        = crit;
    resultsprony.dc1{k}         = dc1;
    resultsprony.dc2{k}         = dc2;
end

end