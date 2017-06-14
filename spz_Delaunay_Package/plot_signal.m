function plot_signal(this_Sig, this_Time_Limit, this_Signal_Limit, Plot_Title, Sample_Rate)
%  plot_signal(this_Sig, this_Time_Limit, this_Signal_Limit, Plot_Title, Sample_Rate)
%   Plots a temporal Signal;
%
% P. Flandrin & Ph. Depalle
% 2015, July 4rd
%
% inputs
%   this_Sig             = temporal signal
%   this_Time_Limit      = [T_Min T_Max] in samples. 
%   this_Signal_Limit    = [A_Min A_Max] in amplitude. 
%   Plot_Title           = Title of the plot
%   Sample_Rate          = Sample_Rate  = 0 -> Time axis in samples
%                          Sample_Rate ~= 0 -> Time axis in seconds,
%                          scaled by Sample_Rate

T_Min = this_Time_Limit(1);
T_Max = this_Time_Limit(2);

this_Time = (T_Min:T_Max);

if(Sample_Rate == 0)
    this_Time_Label = 'Time (in samples)';
else
    this_Time_Label = 'Time (s)';
    this_Time       = this_Time/Sample_Rate;
end

plot(this_Time, this_Sig(T_Min:T_Max));
xlim([this_Time(1) this_Time(end)]);
ylim([this_Signal_Limit(1) this_Signal_Limit(2)]);

xlabel(this_Time_Label);
ylabel('Amplitude');
title(Plot_Title, 'Interpreter', 'none');

drawnow
end