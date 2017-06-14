function [Diff_Sig, SNR_Out] = spz_delaunay_compare(Sig_1, this_Name, Border_Length, SNR_In, Sig_2, Fig_Num, SR_Plot)
% [Diff_Sig, SNR_Out] = spz_delaunay_compare(Sig_1, this_Name, Border_Length, SNR_In, Sig_2, Fig_Num)
%   Computes difference between original signal and reconstructed one.
%   Returns difference signal and SNR.
%
% P. Flandrin & Ph. Depalle
% 2015, June 30th
%
% inputs
%   Sig_1         = Original signal
%   this_Name     = Name of the input signal (for plot's title)
%   Border_Length = Removes first and last Border_Length samples from the comparison
%   SNR_In        = Input noise level in dB
%   Sig_2         = Reconstructed Signal
%   Fig_Num       = Figure index for difference signal
%   SR_Plot       = SR_Plot  = 0 -> Time axis in samples
%                  SR_Plot != 0 -> Time axis in seconds
%
% output
%   Diff_Sig      = Difference between reconstructed signal and original one
%   SNR_Out       = Signal to Noise Ratio between reconstructed signal and original one
%
% Calls:
%   plot_signal
%
   
    fprintf('\nDifference Computation.\n');
    
    Signal_Length = length(Sig_1);
    %Compare_Length = 1:Signal_Length - 0*Nx; % Default domain for comparison
    Compare_Length = ceil(Border_Length):Signal_Length-ceil(Border_Length); % Remove side-effects
    %Compare_Length = 1:Signal_Length;
    
    Diff_Sig = Sig_1(Compare_Length) - Sig_2(Compare_Length);
    max_Diff_Sig = max(abs(Diff_Sig));

    figure(Fig_Num);
    plot_signal(Diff_Sig, [1 (Compare_Length(end) - Compare_Length(1))], [-max_Diff_Sig max_Diff_Sig], [this_Name ' - Difference signal'], SR_Plot);

    SNR_Out = -10*log10(sum(Diff_Sig.^2)/sum(Sig_1(Compare_Length).^2));
    fprintf('\nSNR In: %f \tSNR Out: %f\n', SNR_In, SNR_Out);
end