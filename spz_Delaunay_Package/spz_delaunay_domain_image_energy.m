function Energy_All = spz_delaunay_domain_image_energy(this_Name, this_Tfr, triZ, xZ, yZ, these_Domains, Fig_Num, SR_Plot)
%  Energy_All = spz_delaunay_domain_image_energy(this_Name, this_Tfr, triZ, xZ, yZ, these_Domains, Fig_Num)
%   Returns Energy array for all Domains ; Plots energy as a function of domains and TF energetic regions.
%
% P. Flandrin & Ph. Depalle
% 2015, July 1st
%
% inputs
%   this_Name     = Name of the input signal (for the plot title)
%   this_tfr      = STFT (output from spz_delaunay_dom4a)
%   triZ          = Vector of selected triangles
%   xZ            = Vector of X vertices' coordinates
%   yZ            = Vector of Y vertices' coordinates
%   these_Domains = Cell array of extracted domains
%   Fig_Num       = Figure for selected domains' energy plots.
%   SR_Plot       = SR_Plot  = 0 -> Time axis in samples
%                   SR_Plot ~= 0 -> Time axis in seconds
%
% output
%   Energy_All    = Array Energy of Domains.
%
% Calls:
%   Draw_Sub_Domain
%   plot_tf_delaunay_triangulation
%   plot_tf_image

N_Domains   = length(these_Domains);
[F_Size, T_Size] = size(this_Tfr);

Image_All   = zeros( (F_Size/2) + 1, T_Size);
Energy_All  = zeros(N_Domains, 1);
All_Borders = zeros(N_Domains, 4); % Each row is x_min, x_max, y_min, y_max

fprintf('\n\nCompute Domain Energy.');
for this_Domain = 1:N_Domains
    fprintf('\n\tCompute Domain Energy: %d out of %d.', this_Domain, N_Domains);
    [this_Image, this_Borders] = Draw_Sub_Domain(triZ, xZ, yZ,...
                                                 these_Domains{this_Domain});
    
    All_Borders(this_Domain, :) = this_Borders;
    x_min = this_Borders(1);
    x_max = this_Borders(2);
    y_min = this_Borders(3);
    y_max = this_Borders(4);
    
    this_Sub_Tfr = this_Image .* this_Tfr(x_min:x_max, y_min:y_max);
    this_Energy  = sum(sum( abs(this_Sub_Tfr).^2 ));
    Energy_All(this_Domain) = this_Energy;
    
    this_Image = 10*log10(this_Energy) * this_Image;
    Image_All(x_min:x_max, y_min:y_max) = Image_All(x_min:x_max, y_min:y_max) + this_Image;
end

Energy_Max = max(10*log10(Energy_All));

figure(Fig_Num);
colormap('jet');

subplot(211);
plot_energy(10*log10(Energy_All), [1 N_Domains], [(Energy_Max - 100) (Energy_Max)], [this_Name ' - Energy as a function of the domain index']);

subplot(212);
plot_tf_image(Image_All, [1 T_Size], [1 ((F_Size/2) + 1)], [this_Name ' - Energy of selected Delaunay domains (in dB)'], SR_Plot);
colorbar

fprintf('\n');
end
