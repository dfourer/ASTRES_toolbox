clear all
close all


%  Reassigning and synchrosqueezing the Stockwell Transform
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Mar. 2016

gamma_K = 10^(-4);

selected_signal = 1;  %%choose the signal to process: can be 1 or 2
figs_folder = 'figs';

if selected_signal == 1
    load fig1-observed-Livingston.txt
    t1=fig1_observed_Livingston(:,1);
    sig1=fig1_observed_Livingston(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    %figs_folder = 'figs_signal1';
elseif selected_signal == 2

    load fig1-observed-Hanford.txt
    t1=fig1_observed_Hanford(:,1);
    sig1=fig1_observed_Hanford(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    %figs_folder = 'figs_signal2';
else
    error('Unknown signal')
end

if ~exist(figs_folder, 'dir')
 mkdir(figs_folder);
end

NbPoints=length(sig1);   %length(sig1)
alpha=0.5;             %% TFR compression factor (modify display contrast)

%% TF analysis
w0T = 8;
w0  = 2 * pi * 100;   %% used for wavelet

T = w0T/w0;           %% used for gabor transform
L = T/Ts;


M=1000; %4000; %5000     %% number of frequency bins

mi=floor(20*M/Fe1); mf=ceil(520*M/Fe1);   %% frequency range [20Hz - 520Hz]
s = sig1;
s = s-mean(s);  % center signal

t = t1;
m = m_axis(M);
nfreqs = m/M;

m_range = mi:mf;



%% compute the synchrosqueezed S-transform
[tfr, stfr, lost] = tfrsst(s, M, w0T, gamma_K);

%[tfr, stfr, lost] = tfrvsst(s, M, w0T, gamma_K);  %%vertical synchrosqueezed ST







% 
% figure
% imagesc(t, nfreqs(m_range), (abs(stfr(m_range,:)).^2).^alpha);
% set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap)); 
% xlabel('time (s)', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
% ylabel('normalized frequency', 'FontSize', 16)
% title(sprintf('Vertical synchrosqueezed S-transform, w0T=%.2f', w0T), 'FontSize', 14);
% 
% 
% %% TF filtering 
% %hist(abs(reshape(stfr,1,size(stfr,1) * size(stfr,2))),100)
%    
% %% test1: global hard thresholding
% T = 0.2;
% mask1 = zeros(size(stfr));
% mask1(abs(stfr) > T) = 1;
% 
% 
% 
% %%manual adjustments of TF-mask
% 
% mask1(13:end, 1:2000) = 0;      %%% all frequencies above 50Hz are set to 0 until sample index n=1800
% mask1(25:end, 2000:2500) = 0;   %%% all frequencies above 65Hz are set to 0 for sample index in [1800:2700]
% mask1(1:35, 2900:end) = 0;  
% %mask1(75:end, 2900:end) = 0;
% %mask1(1:12, 2900:end) = 0;
% mask1(8:13, 1880:2000) = 1;
% 
% stfr_filtered = stfr .* mask1;
% 
% imagesc((abs(stfr_filtered(m_range,:))))
% set(gca,'YDir','normal');
% 
% 
% %% test2: local hard thresholding based on the median
% %a = 2;   %% threshold factor
% %stfr_filtered = stfr;
% %for i = 1:size(stfr_filtered,2)
% %     
% %  T = a * mean(abs(stfr_filtered(:,i)));
% %  I = find(abs(stfr_filtered(:,i)) < T);
% %  stfr_filtered(I,i ) = 0;
% %end
%     
% %display the difference before TF-filtering process
% alpha = 0.3;
% figure(2)
% subplot(121)
% imagesc(abs(stfr(m_range,:)).^(2*alpha))
% set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap)); 
% set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
% title('before filtering')
% subplot(122)
% imagesc(abs(stfr_filtered(m_range,:)).^(2*alpha))
% set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap));
% set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
% title('after filtering')
% %subplot(133)
% %imagesc(abs(stfr(m_range,:)-stfr_filtered(m_range,:)).^(2*alpha))
% %set(gca,'YDir','normal');colormap gray;colormap(flipud(colormap));
% %set(gca,'YTick',[]);set(gca,'XTick',[]);xlabel('time');ylabel('frequency')
% %title('difference')
%     
% 
% figure(3);
% imagesc(t, nfreqs(m_range), (abs(stfr_filtered(m_range,:)).^2).^alpha);
% set(gca,'YDir','normal')
% colormap gray;colormap(flipud(colormap)); 
% xlabel('time (s)', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
% ylabel('normalized frequency', 'FontSize', 16)
% title(sprintf('Filtered Vertical Synchrosqueezed S-Transform w0T=%.2f', w0T), 'FontSize', 14);
% %    saveas(gcf, sprintf('%s/rsst_%d.eps', figs_folder, index), 'epsc');
% 
% %n_range_rec = n_range-n0;
% s_hat = rectfrsst(stfr_filtered(m_range,:), w0T, m_range);
% s_rec = 2*real(s_hat); %% assume the signal is real
% 
% figure
% plot(t, s, 'g-.');  %t(n_range_rec), 
% hold on
% plot(t, s_rec, 'k'); %t(n_range_rec), 
% xlabel('time (s)', 'FontSize', 16)
% legend('reference signal', 'reconstructed signal');
% title('Resulting signal')
