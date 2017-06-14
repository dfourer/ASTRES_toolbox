function X_Out = spz_delaunay_reconstruct(this_Signal, this_Name, this_tfr, triZ, xZ, yZ, this_Domains, Dom_Select, Mask_Flag, Fig_Num, SR_Plot)
% X_Out = spz_delaunay_reconstruct(this_Signal, this_Name, this_tfr, triZ, xZ, yZ, this_Domains, Dom_Select, Mask_Flag, Fig_Num)
%   Returns signal reconstructed from the selected domains
%
% P. Flandrin & Ph. Depalle
% 2015, June 30th
%
% inputs
%   this_Signal  = Input signal
%   this_Name    = Name of the input signal (for the plot title)
%   this_tfr     = STFT (output from spz_delaunay_dom4a)
%   triZ         = Vector of selected triangles
%   xZ           = Vector of X vertices' coordinates
%   yZ           = Vector of Y vertices' coordinates
%   this_Domains = Cell array of extracted domains
%   Dom_Select   = Index array of selected domains (from this_Domains)
%   Mask_Flag    = true - use selected domains; false - use transparent mask
%   Fig_Num      = Figure index for mask and reconstructed signal
%                  Fig_Num+1 for original and reconstructed signal
%   SR_Plot      = SR_Plot  = 0 -> Time axis in samples
%                  SR_Plot != 0 -> Time axis in seconds
% output
%   X_Out        = reconstructed signal from selected domains
%
% Calls:
%   Draw_Sub_Domain
% 	get_Plural
% 	plot_signal
% 	plot_tf_image
% 	roundgauss
% 	spz_delaunay_dom4c2

fprintf('\n\nReconstruction.\n');

[Nx,Signal_Length]  = size(this_tfr) ;
prec                = 10^(-6) ;
[H, ~]              = roundgauss(Nx, prec);
max_window          = max(H);
max_this_Signal     = max(abs(this_Signal));

Nx              = Nx/2;
NAx             = Nx + 1;
NAy             = Signal_Length;

N_Domains       = length(this_Domains);
domcol          = N_Domains + 1 - Dom_Select;
in_the_Plural   = get_Plural(length(domcol));


if(Mask_Flag)
    mask = zeros(NAx, NAy) ;
    Title_Mask   = [this_Name ' - Mask from selected domain'...
                           in_the_Plural];
    Title_Signal = [this_Name ' - Reconstructed component' in_the_Plural ' from mask'];
                       
    for this_color = 1:length(domcol)
        this_dom = domcol(this_color);
        [this_Mask, this_Borders] = Draw_Sub_Domain(triZ, xZ, yZ, this_Domains{this_dom});
        
        x_min = this_Borders(1);
        x_max = this_Borders(2);
        y_min = this_Borders(3);
        y_max = this_Borders(4);
        mask(x_min:x_max, y_min:y_max) = mask(x_min:x_max, y_min:y_max) + this_Mask;
    end
else
    mask         = ones(NAx, NAy);
    Title_Mask   = [this_Name ' - Transparent mask'];
    Title_Signal = [this_Name ' - Reconstructed from transparent mask'];
end

X_Out   = spz_delaunay_dom4c2(this_tfr, mask) / (max_window * Nx);
X_Out   = real(X_Out);

figure(Fig_Num)
subplot(211)
plot_tf_image(mask, [1 Signal_Length], [1 NAx], Title_Mask, SR_Plot);

subplot(212)
plot_signal(X_Out, [1 Signal_Length], [-max_this_Signal max_this_Signal], Title_Signal, SR_Plot);


figure(Fig_Num+1)
subplot(211);
plot_signal(this_Signal, [1 Signal_Length], [-max_this_Signal max_this_Signal], [this_Name ' - Original signal'], SR_Plot);

subplot(212);
plot_signal(X_Out, [1 Signal_Length], [-max_this_Signal max_this_Signal], [this_Name ' - Reconstructed signal'], SR_Plot);
end
