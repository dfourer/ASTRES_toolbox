%%   Script used to generate all figures
%
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 24-07-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassigning and Synchrosqueezing of the Stockwell Transform. IEEE. Tran. on Signal Processing. 2015]
clear all
close all

NbImg  = 5;
NbImg2 = 5; %% for LMR

signal = 1; %5;  %% change value in [1,2,3,4]


mu_range  = [0.01 0.05 0.2 0.5 1 2 5];      %mu_range = logspace(log10(mu_range(1)), log10(mu_range(2)), NbImg2);
SNR_range = [45 35 25 15 5];
%w0T_range = [2*pi*2 2*pi*5 2*pi*10];        %w0T_range = linspace(2*pi*1, 2*pi*15, NbImg); % (2*pi*10)          %% range for omega_0 T
w0T_range = 2*pi*[1.5 2 3.5];


alpha = 0.3;

Fs = 1000;          %% sampling frequency
Ts = 1/Fs;          %% sampling period
%t  = 0:Ts:(1-Ts);   %% time axis
M  = 500;           %% number of frequency bins
M2 = round(M/2);

load_signal
n  = 1:length(s);   %% time sample
N  = length(n);     %% signal length

s0 = s-mean(s); %% 0-mean signal
s = s0;

%%% prepare folders to store images
chemin0 = 'figs_TSP';

cheminr = sprintf('%s/%d', chemin0, signal);
chemin{1} = sprintf('%s/st/',  cheminr);
index1 = 1;
chemin{2} = sprintf('%s/rst/', cheminr);
index2 = 1;
chemin{3} = sprintf('%s/sst', cheminr);
index3 = 1;
chemin{4} = sprintf('%s/vsst', cheminr);
index4 = 1;

for i = 1:length(chemin)
 if ~exist(chemin{i}, 'dir')
   mkdir(chemin{i});
 end
end

gamma_K = 10^(-4);
f      = ((1:M2)-1)/M * Fs;
nfreqs = f / Fs;

fp = fopen(sprintf('%s/results.txt',cheminr), 'w+b');

for id = 1:length(SNR_range)
    
 %% add noise
 rsb_target = SNR_range(id);
 s = sigmerge(s0, randn(size(s0)), rsb_target);

 rsb = SNR(s0, s);
  if abs(rsb - rsb_target) > 10
    error('Invalid SNR')
  end
 
 tmp = sprintf('+++++++++++++  SNR= %0.02f +++++++++++++\n', rsb);
 fprintf(1, tmp); fprintf(fp, tmp);
   
 for idx1 = 1:length(w0T_range)
  w0T = w0T_range(idx1);
  f0T = w0T/(2*pi);
  
  tmp = sprintf('++ W0T=%0.02f, f0T=%0.02f ++\n', w0T, f0T);
  fprintf(1, tmp); fprintf(fp, tmp);
 
  %% 1 Compute the classical S-transform and display Stockwellogram
  %[tfr] = tfrst(s, M, w0T, gamma_K);
  [tfr, rtfr] = tfrrst(s, M, w0T, gamma_K);
  figure(1)
  imagesc(n, nfreqs, (abs(tfr(1:M2,:)).^2).^alpha);
  %set(gca,'Ytick',0:0.05:0.5)
  %set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Stockwellogram SNR=%0.02f, f_0T=%0.02f', rsb_target, f0T), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/st_%d.eps', chemin{1}, index1), 'epsc');
  index1 = index1+1;
 
  %% 1a) Reconstruction classique
  [s_hat] = rectfrst(tfr);
  tmp = sprintf('S-transform reconstruction (1st formula), RQF=%0.02f, f_0 T=%0.02f\n', RQF(s,s_hat'), f0T);
  fprintf(1, tmp); fprintf(fp, tmp);
 
  %% 1b) Simplified reconstruction formula (OK)
  [s_hat2] = rectfrst2(tfr, w0T);
  tmp = sprintf('S-transform reconstruction (2nd formula), RQF=%0.02f, f_0 T=%0.02f\n', RQF(s,s_hat2'), f0T);
  fprintf(1, tmp); fprintf(fp, tmp);
 
  
  %% 2 Compute the Reassigned Stockwellogram
  %[tfr, rtfr] = tfrrst(s, M, w0T, gamma_K);
  figure(2)
  imagesc(n, nfreqs, (abs(rtfr(1:M2,:)).^2).^alpha);
  %set(gca,'Ytick',0:0.05:0.5)
  %set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Reassigned Stockwellogram SNR=%0.02f, f_0T=%0.02f', rsb_target, f0T), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/rst_%d.eps', chemin{2}, index2), 'epsc');
  index2 = index2+1;
 
 %% 3 Compute the Synchrosqueezed S-transform
 [~,stfr] = tfrsst(s, M, w0T, gamma_K);
 figure(3)
 imagesc(n, nfreqs, (abs(stfr(1:M2,:)).^2).^alpha);
 %set(gca,'Ytick',0:0.05:0.5)
 %set(gca,'Xtick',0:50:N)
 ylim([0 0.499])
 set(gca,'YDir','normal')
 xlabel('time samples', 'FontSize', 16)
 ylabel('normalized frequency', 'FontSize', 16)
 title(sprintf('Synchrosqueezed ST SNR=%0.02f, f_0T=%0.02f', rsb_target, f0T), 'FontSize', 14);
 colormap gray;
 cmap = colormap;
 cmap = flipud(cmap);
 colormap(cmap);
 saveas(gcf, sprintf('%s/sst_%d.eps', chemin{3}, index3), 'epsc');
 index3 = index3+1;
 
 %% 3b) Reconstruct the Synchrosqueezed S-transform
 %s_hat3 = real(rectfrsst(stfr, w0T)); fprintf(1,'Synchrosqueezed ST reconstruction, SNR=%0.02f, f0 T=%0.02f\n', SNR(s,s_hat3'), f0T);
 
 %s_hat3 = 2*real(rectfrsst(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
 s_hat3 = rectfrsst(stfr, w0T);
 
 tmp = sprintf('Synchrosqueezed ST reconstruction, SNR=%0.02f, f0 T=%0.02f\n', RQF(s,s_hat3'), f0T);
 fprintf(1, tmp); fprintf(fp, tmp);
 
  %% 4) Compute the Levenberg-Marquardt Reassigned Stockwellogram
 for idx2 = 1:length(mu_range)
  mu = mu_range(idx2);
  [~, rtfr, lost] = tfrlmrst(s, mu, M, w0T, gamma_K);
  figure(2)
  imagesc(n, nfreqs, (abs(rtfr(1:M2,:)).^2).^alpha);
  %set(gca,'Ytick',0:0.05:0.5)
  %set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('LM-Reassigned Stockwellogram SNR=%0.02f, f_0T=%0.02f \\mu=%0.02f', rsb_target, f0T, mu), 'FontSize', 14);
  %colormap gray;
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/rst_%d.eps', chemin{2}, index2), 'epsc');
  index2 = index2+1;
 end
 

  %% 5) Compute LM synchrosqueezed ST
 for idx2 = 1:length(mu_range)
  mu = mu_range(idx2);
  [~, stfr, lost] = tfrlmsst(s, mu, M, w0T, gamma_K);
  figure(3)
  imagesc(n, nfreqs, (abs(stfr(1:M2,:)).^2).^alpha);
  %set(gca,'Ytick',0:0.05:0.5)
  %set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('LM Synchrosqueezed ST SNR=%0.02f dB, f_0T=%0.02f \\mu=%0.02f', rsb_target, f0T, mu), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/sst_%d.eps', chemin{3}, index3), 'epsc');
  index3 = index3+1;
  
  s_hat3 = rectfrsst(stfr, w0T);
  tmp = sprintf('Synchrosqueezed ST reconstruction, RQF=%0.02f dB, mu=%0.02f, f0T=%0.02f\n', RQF(s,s_hat3'), mu, f0T);
  fprintf(1, tmp); fprintf(fp, tmp);
 end
 
 %% 6) Compute the Vertical second-order synchsqueezed ST
 [~,stfr] = tfrvsst(s, M, w0T, gamma_K);
 figure(4)
 imagesc(n, nfreqs, (abs(stfr(1:M2,:)).^2).^alpha);
 %set(gca,'Ytick',0:0.05:0.5)
 %set(gca,'Xtick',0:50:N)
 ylim([0 0.499])
 set(gca,'YDir','normal')
 xlabel('time samples', 'FontSize', 16)
 ylabel('normalized frequency', 'FontSize', 16)
 title(sprintf('Vertical Synchrosqueezed ST SNR=%0.02f, f_0T=%0.02f', rsb_target, f0T), 'FontSize', 14);
 colormap gray;
 cmap = colormap;
 cmap = flipud(cmap);
 colormap(cmap);
 saveas(gcf, sprintf('%s/vsst_%d.eps', chemin{4}, index4), 'epsc');
 index4 = index4+1;
 
 end    %% w0T
end  %% SNR

fclose all;
for i = 1:length(chemin)
 eps2pdf(chemin{i});
end