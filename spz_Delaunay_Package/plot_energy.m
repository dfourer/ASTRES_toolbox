function plot_energy(this_Energy, this_Index_Limit, this_Energy_Limit, this_Energy_Title)

if(length(this_Energy) ~= 1)
    Index_Min = this_Index_Limit(1);
    Index_Max = this_Index_Limit(2);
    plot(Index_Min:Index_Max, this_Energy(Index_Min:Index_Max));
    xlim([Index_Min Index_Max]);
    ylim([this_Energy_Limit(1) this_Energy_Limit(2)]);
    title(this_Energy_Title, 'Interpreter', 'none');
else
    Energy_Max = this_Energy_Limit(2);
    plot(this_Energy, '+', 'MarkerSize', 16, 'MarkerEdgeColor', 'b');
    xlim([0.999 1.001]);
    ylim([(Energy_Max - 50) (Energy_Max + 50)]);
    title([this_Energy_Title ' (only ONE value)'], 'Interpreter', 'none');
end

xlabel('Domain index');
ylabel('Energy (dB)');

drawnow
end