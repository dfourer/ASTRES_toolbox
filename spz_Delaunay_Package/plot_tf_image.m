function plot_tf_image(this_Image, this_Time_Limit, this_Frequency_Limit, Plot_Title, Sample_Rate)
%  plot_tf_image(this_Image, this_Time_Limit, this_Frequency_Limit, Plot_Title, Sample_Rate)
%   Plots a Time-Frequency Image;
%
% P. Flandrin & Ph. Depalle
% 2015, July 4rd
%
% inputs
%   this_Image           = Time Frequency Image; 1st dimension is frequency, 2nd domension time
%   this_Time_Limit      = [T_Min T_Max] in samples. 
%   this_Frequency_Limit = [F_Min F_Max] in frequency bins. 
%   Plot_Title           = Title of the plot
%   Sample_Rate          = Sample_Rate  = 0 -> Time axis in samples
%                          Sample_Rate ~= 0 -> Time axis in seconds,
%                          scaled by Sample_Rate

[F_Size, ~] = size(this_Image); % Assumes F_Size is 1+ Nfft/2
N_FFT = 2*(F_Size - 1);

T_Min = this_Time_Limit(1);
T_Max = this_Time_Limit(2);

F_Min = this_Frequency_Limit(1);
F_Max = this_Frequency_Limit(2);

this_Time       = (T_Min:T_Max);
this_Frequency  = (F_Min:F_Max)/N_FFT;

if(Sample_Rate == 0)
    this_Time_Label      = 'Time (in samples)';
    this_Frequency_Label = 'Normalized frequency';
else
    this_Time_Label      = 'Time (s)';
    this_Frequency_Label = 'Frequency (Hz)';
    this_Time       = this_Time      / Sample_Rate;
    this_Frequency  = this_Frequency * Sample_Rate;
end

imagesc(this_Time, this_Frequency, this_Image(F_Min:F_Max, T_Min:T_Max));
set(gca,'YDir','normal');
xlim([this_Time(1)      this_Time(end)]);
ylim([this_Frequency(1) this_Frequency(end)]);

xlabel(this_Time_Label);
ylabel(this_Frequency_Label);
title(Plot_Title, 'Interpreter', 'none');

drawnow
end