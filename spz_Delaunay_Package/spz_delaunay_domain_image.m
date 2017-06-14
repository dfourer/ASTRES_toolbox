function Image_All = spz_delaunay_domain_image(this_Name, triZ, xZ, yZ, these_Domains, X_Size, Y_Size, Fig_Num, SR_Plot)
%  Image_All = spz_delaunay_domain_image(this_Name, triZ, xZ, yZ, these_Domains, X_Size, Y_Size, Fig_Num)
%  Returns Domain image; Plots selected domains' triangles and image.
%
%  P. Flandrin & Ph. Depalle
%  2015, July 1st
%
% inputs
%   this_Name     = Name of the input signal (for the plot title)
%   triZ          = Vector of selected triangles
%   xZ            = Vector of X vertices' coordinates
%   yZ            = Vector of Y vertices' coordinates
%   these_Domains = Cell array of extracted domains
%   X_Size        = Number of spectral bins (Nx+1)
%   Y_Size        = Number of time samples (Signal_Length)
%   Fig_Num       = Figure for selected domains' triangles and image plots.
%   SR_Plot       = SR_Plot  = 0 -> Time axis in samples
%                   SR_Plot ~= 0 -> Time axis in seconds
%                
% output
%   Image_All     = domain image.
%
% Calls:
%   Draw_Sub_Domain
%   plot_tf_delaunay_triangulation
%   plot_tf_image

N_Domains   = length(these_Domains);

Image_All   = zeros(X_Size,Y_Size);
All_Borders = zeros(N_Domains, 4); % Each row is x_min, x_max, y_min, y_max

fprintf('\n\nCompute Domain Image.');
for this_domain = 1:N_Domains
    fprintf('\n\tCompute Domain Image: %d out of %d.', this_domain, N_Domains);
    [this_Image, this_Borders] = Draw_Sub_Domain(triZ, xZ, yZ,...
        these_Domains{this_domain});
    
    All_Borders(this_domain, :) = this_Borders;
    x_min = this_Borders(1);
    x_max = this_Borders(2);
    y_min = this_Borders(3);
    y_max = this_Borders(4);

    this_Image = (N_Domains + 1 - this_domain) * this_Image;
    Image_All(x_min:x_max, y_min:y_max) = Image_All(x_min:x_max, y_min:y_max) + this_Image;
end

figure(Fig_Num);
clf(Fig_Num);
subplot(211);

[y_scale, x_scale] = plot_tf_delaunay_triangulation(triZ, xZ, yZ, X_Size, [1 Y_Size], [1 X_Size], [this_Name ' - Delaunay triangulation (selected triangles - no pruning)'], SR_Plot);

xZ = x_scale * xZ;
yZ = y_scale * yZ;

this_colormap = colormap('jet');
hold on
for this_domain = 1:N_Domains    
    x_min = x_scale * All_Borders(this_domain, 1);
    x_max = x_scale * All_Borders(this_domain, 2);
    y_min = y_scale * All_Borders(this_domain, 3);
    y_max = y_scale * All_Borders(this_domain, 4);
    
    this_Color = this_colormap(ceil(64*(N_Domains + 1 - this_domain)/N_Domains),:);
    rectangle('Position', [y_min x_min (y_max - y_min) (x_max - x_min)],...
              'EdgeColor', this_Color, 'LineWidth', 3);
end
hold off
drawnow

subplot(212);
plot_tf_image(Image_All, [1 Y_Size], [1 X_Size], [this_Name ' - Selected Delaunay domains'], SR_Plot);
fprintf('\n');
end
