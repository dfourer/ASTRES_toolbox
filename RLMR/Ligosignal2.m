% Recursive versions of the Levenberg-Marquardt reassigned spectrogram 
% and of the synchrosqueezed STFT short demo
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Feb. 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

clear all
close all


selected_signal = 1;  %%choose the signal to process: can be 1 or 2

if selected_signal == 1
    load fig1-observed-Livingston.txt
    t1=fig1_observed_Livingston(:,1);
    sig1=fig1_observed_Livingston(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal1';
elseif selected_signal == 2

    load fig1-observed-Hanford.txt
    t1=fig1_observed_Hanford(:,1);
    sig1=fig1_observed_Hanford(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal2';
else
    error('Unknown signal')
end

if ~exist(figs_folder, 'dir')
 mkdir(figs_folder);
end
NbPoints=length(sig1);   %length(sig1)
alpha=0.5;             %% TFR compression factor (modify display contrast)

%% TF analysis
k = 4;      %% recursive filter order
L = 60;     %% analysis window length

M=5000;     %% number of frequency bins
mi=floor(20*M/Fe1); mf=ceil(520*M/Fe1);   %% frequency range [20Hz - 520Hz]

n0     = (k-1) * L;
N_0pad = n0 + 20;                           % N_0pad should be greater than n0
AddedZeros = zeros(N_0pad, 1);
n_range  = N_0pad+(1:length(sig1));
s = real([AddedZeros;sig1;AddedZeros]);     
N = length(sig1); %NbPoints+N_0pad;         %% signal length 
t = 1:length(s)-N_0pad;                     %% time index without the first added zeros


%% compute the recursive synchrosqueezed transform
[~, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0);
figure(1)
imagesc(t, nfreqs, (abs(stfr(:,n_range)).^2).^alpha);
set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Recursive synchrosqueezed STFT k=%d, L=%2.2f ', k, L), 'FontSize', 14);
    
    
%% TF filtering 
%hist(abs(reshape(stfr,1,size(stfr,1) * size(stfr,2))),100)
   
%% test1: global hard thresholding
T = 0.2;
mask1 = zeros(size(stfr));
mask1(abs(stfr) > T) = 1;

%%manual adjustments of TF-mask
mask1(15:end, 1:1800) = 0;      %%% all frequencies above 50Hz are set to 0 until sample index n=1800
mask1(20:end, 1800:2700) = 0;   %%% all frequencies above 65Hz are set to 0 for sample index in [1800:2700]
mask1(30:end, 2700:2900) = 0;  
mask1(75:end, 2900:end) = 0;
mask1(1:12, 2900:end) = 0;

stfr_filtered = stfr .* mask1;
  

%% test2: local hard thresholding based on the median
%a = 2;   %% threshold factor
%stfr_filtered = stfr;
%for i = 1:size(stfr_filtered,2)
%     
%  T = a * mean(abs(stfr_filtered(:,i)));
%  I = find(abs(stfr_filtered(:,i)) < T);
%  stfr_filtered(I,i ) = 0;
%end
    
%display the difference before TF-filtering process
alpha = 0.3;
figure(2)
subplot(131)
imagesc(abs(stfr).^(2*alpha))
set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap)); 
set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
title('before filtering')
subplot(132)
imagesc(abs(stfr_filtered).^(2*alpha))
set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap));
set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
title('after filtering')
subplot(133)
imagesc(abs(stfr-stfr_filtered).^(2*alpha))
set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap));
set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
title('difference')
    

figure(3);
imagesc(t, nfreqs, (abs(stfr_filtered(:,n_range)).^2).^alpha);
set(gca,'YDir','normal')
colormap gray;colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Filtered Recursive synchrosqueezed STFT k=%d, L=%2.2f ', k, L), 'FontSize', 14);
%    saveas(gcf, sprintf('%s/rsst_%d.eps', figs_folder, index), 'epsc');

n_range_rec = n_range;
%n_range_rec = n_range-n0;
s_hat = sstft_rec(stfr_filtered, k, L, M, n0);
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
    
figure
plot(t(n_range_rec)*Te1, s(n_range_rec), 'g-.');  %t(n_range_rec), 
hold on
plot(t(n_range_rec)*Te1, s_rec, 'k'); %t(n_range_rec), 
xlabel('time (s)', 'FontSize', 16)
legend('reference signal', 'reconstructed signal');
title('Resulting signal')
