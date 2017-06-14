function [tfr, triZ_Final, xZ_Final, yZ_Final] = spz_delaunay_dom4a(x, Nx, Example_Name, Fig_Num, SR_Plot)
%
% P. Flandrin, May 2015
%
% input
%   x = signal
%   Nx = Block Size (has to be a power of 2, and preferably <= 1024)
%   fignum = figure number for Delaunay Triangulation
%   SR_Plot       = SR_Plot  = 0 -> Time axis in samples
%                   SR_Plot ~= 0 -> Time axis in seconds
%
% output
%   tfr                   STFT
%   triZ_Final            Set of triangles (border effects removed)
%   xZ_Final, yZ_Final    Coordinates of triangles vertices
%
% needs
%   Time-Frequency ToolBox from tftb.nongnu.org
% calls
%   roundgauss.m
%   tfr_stft.m
%   extr2minth.m
%   Catenate_Tri.m
%   Select_Triangles.m
%   plot_tf_delaunay_triangulation
%	plot_tf_image
 


%%  
   
Nfft = 2*Nx ;%2*Nx ; % number of FFT bins

prec = 10^(-6) ;
[H, L] = roundgauss(Nfft, prec);

t_border = 2*ceil(L);
ttd = t_border + 1;
ttf = Nx - t_border;
ffd = ceil(L/2) + 1;
fff = Nfft/2 - ceil(L/2);

% Offset_Freq = 2*Nx/64;
% Offset_Time = Nx/2/c;
Offset_Freq = ceil((Nfft/2)/L);

 
Signal_Length   = length(x);
Begin_Time      = 0;    % Should be zero I think
Block_Duration  = Nx;
Delta_T         = ttf - ttd;
Step_Incr       = Block_Duration - 2*t_border;

% % June 24th: Modify x's length so that all signal is being anlysed
% 
% Begin_Time_Max  = (Signal_Length - Nx/64 - Block_Duration);
% x               = [x; zeros(mod(Begin_Time_Max, Step_Incr), 1)];
% size(x)
% Signal_Length   = length(x);


% spectro and STFT
% 
% [tfr,rtfr,hat] = tfrrsp_sq(x,t,Nfft,H,1) ;
t   = 1:Signal_Length ;
tfr = tfr_stft(x, t, Nfft, H, 1);
 
S_0 = abs(tfr);
Lx = Nfft/L ;
Ly = L ;

% all - Build S

Total_Freq = Nfft/2 + 2*Offset_Freq;
S = zeros(Total_Freq, Signal_Length);
S(1+Offset_Freq:Total_Freq,:) = S_0(1:Total_Freq - Offset_Freq,:);
S(1:Offset_Freq, :) =  S_0(Nfft-Offset_Freq+1:Nfft, :);


triZ_Final      = [];
xZ_Final        = [];
yZ_Final        = [];
last_min_index  = 1;  % This one is just to speed up triangle catenation

%---- Plot STFT of the signal - BEGIN
figure(Fig_Num)  % Start here to allow drawing rectangles within the while loop below
clf(Fig_Num)     % To prevent superimposition between two experiments
subplot(211)

S_Max = 20*log10(max(max(S_0)));
S_Min = 20*log10(min(min(S_0)));
Color_Min =  1 - (80/(S_Max - S_Min));
if(Color_Min < 0)
    Color_Min = 0;
end

plot_tf_image((20*log10(S_0(1:Nx+1,:))-S_Min)/(S_Max-S_Min), [1 Signal_Length], [1 Nx+1], [Example_Name ' - STFT Magnitude'], SR_Plot);
caxis([Color_Min 1]);
colormap('jet');
drawnow
%---- Plot STFT of the signal - END


%Step = [Nx/64; Block_Duration - 64];
Begin_Time = Nx/64;
first_time = true;
Rectangle_Index = 0;
fprintf('Begin Time:\n');
while(Begin_Time+Block_Duration <= Signal_Length)
    fprintf('\n\tBegin Time: %d End Time: %d', Begin_Time, Begin_Time+Block_Duration);
    Sl = S(:, Begin_Time+1:Begin_Time+Block_Duration);
    [xZl,yZl] = extr2minth(Sl, max(max(abs(Sl)))/10^14); 
    triZl     = delaunay(xZl/Lx,yZl/Ly) ;
    fprintf('\tVertices: %d,\tTriangles: %d, \t Thres.: %e',...
            length(xZl), length(triZl), max(max(abs(Sl)))/10^14);
    % Remove triangle which spans over time borders.

    triZafl = Select_Triangles(triZl, xZl, yZl, ffd,...
                      1 + (Total_Freq - Offset_Freq),ttd, ttf, first_time);
    % Restore proper time/frequency coordinates
    yZl = yZl + Begin_Time;
    xZl = xZl - Offset_Freq;
    
    % Catenat current triangle selection to global one.
    [triZ_Final, xZ_Final, yZ_Final, min_index] =...
        Catenate_Tri(triZ_Final, xZ_Final, yZ_Final, triZafl, xZl, yZl, last_min_index);
    last_min_index = min_index;
    
    % Store Block positions for future rectangle plots
    Rectangle_Index          = Rectangle_Index +1;
    Rect(Rectangle_Index, 1) = ttd+Begin_Time;
    Rect(Rectangle_Index, 2) = 1;
    Rect(Rectangle_Index, 3) = Delta_T;
    Rect(Rectangle_Index, 4) = Nfft/2;

    Begin_Time  = Begin_Time + Block_Duration - 2*t_border; %- Nx/c;
    first_time = false;
end
fprintf('\n');

subplot(212);
[t_scale, f_scale] = plot_tf_delaunay_triangulation(triZ_Final, xZ_Final, yZ_Final, Nx+1, [1 (Begin_Time +4*ceil(L))], [1 (Nx+1)], [Example_Name ' - Delaunay triangulation'], SR_Plot);

% Rectangle plots
hold on
for this_Rectangle_Index = 1:Rectangle_Index
    this_Rectangle = Rect(this_Rectangle_Index, :);
    this_Rectangle([1 3]) = t_scale * this_Rectangle([1 3]);
    this_Rectangle([2 4]) = f_scale * this_Rectangle([2 4]);
    rectangle('Position', this_Rectangle, 'EdgeColor', 'k', 'LineWidth', 3)
end
hold off
drawnow
