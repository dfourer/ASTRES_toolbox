function [this_T_Scale, this_F_Scale] = plot_tf_delaunay_triangulation(this_Triangle, this_xZ, this_yZ, F_Size, this_Time_Limit, this_Frequency_Limit, Plot_Title, Sample_Rate)
%  [this_T_Scale, this_F_Scale] = plot_tf_delaunay_triangulation(this_Triangle, this_xZ, this_yZ, F_Size, this_Time_Limit, this_Frequency_Limit, Plot_Title, Sample_Rate)
%  Plots Delaunay triangulation;
%  Returns Time and Frequency scaling factors for proper scaling of the axes
%
%  P. Flandrin & Ph. Depalle
%  2015, July 6th
%
% inputs
%   this_Triangle        = Vector of selected triangles
%   this_xZ              = Vector of X vertices' coordinates
%   this_yZ              = Vector of Y vertices' coordinates
%   X_Size               = Number of spectral bins (Nx+1)
%   this_Time_Limit      = [T_Min T_Max] in samples. 
%   this_Frequency_Limit = [F_Min F_Max] in frequency bins. 
%   Plot_Title           = Title of the plot
%   Sample_Rate          = Sample_Rate  = 0 -> Time axis in samples
%                          Sample_Rate ~= 0 -> Time axis in seconds,
%                          scaled by Sample_Rate
% output
%   this_T_Scale         = Time scaling factor for the horizontal axis.
%   this_F_Scale         = Frequency scaling factor for the vertical axis.

N_FFT = 2*(F_Size - 1);
 
T_Min = this_Time_Limit(1);
T_Max = this_Time_Limit(2);

F_Min = this_Frequency_Limit(1);
F_Max = this_Frequency_Limit(2);

if(Sample_Rate == 0)
     this_Time_Label      = 'Time (in samples)';
     this_Frequency_Label = 'Normalized frequency';
     this_T_Scale         = 1;
     this_F_Scale         = 1/N_FFT;
else
     this_Time_Label      = 'Time (s)';
     this_Frequency_Label = 'Frequency (Hz)';
     this_T_Scale         = 1/Sample_Rate;
     this_F_Scale         = Sample_Rate/N_FFT;
end

this_Time       = this_T_Scale * (T_Min:T_Max);
this_Frequency  = this_F_Scale * (F_Min:F_Max);
this_yZ         = this_T_Scale * this_yZ;
this_xZ         = this_F_Scale * this_xZ;

triplot(this_Triangle, this_yZ, this_xZ, 'k');

xlim([this_Time(1)      this_Time(end)]);
ylim([this_Frequency(1) this_Frequency(end)]);
 
xlabel(this_Time_Label);
ylabel(this_Frequency_Label);
title(Plot_Title, 'Interpreter', 'none');
drawnow
end