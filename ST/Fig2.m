%%  Script used to generate all figures
%
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 24-07-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassigning and Synchrosqueezing of the Stockwell Transform. IEEE. Tran. on Signal Processing. 2015]
clear all
close all

signal = 1; %5;  %% change value in [1,2,3,4]

%
mu_range  = [0.01 0.05 0.2 0.5 1 2 5];       %mu_range = logspace(log10(mu_range(1)), log10(mu_range(2)), NbImg2);
SNR_range = [45 25 10];                      %[45 35 25 15 5];
%w0T_range = [2*pi*0.8 2*pi*2 2*pi*5];    % 2*pi*10
%w0T_range = [2*pi*2 2*pi*5 2*pi*10];
w0T_range = 2*pi*[1.5 2 3.5];

% 
w0          = 2 * pi * 100;
w0T_default = 7 * pi;

alpha = 0.3; % 0.3;

Fs = 1000;          %% sampling frequency
Ts = 1/Fs;          %% sampling period
%t  = 0:Ts:(1-Ts);   %% time axis
M  = 500;           %% number of frequency bins
M2 = round(M/2);

load_signal
n  = 1:length(s);   %% time sample
N  = length(n);     %% signal length
s0 = s-mean(s); %% 0-mean signal


%%% prepare folders to store images
chemin0 = 'figs_TSP';

cheminr = sprintf('%s/%d', chemin0, signal);
chemin{1} = sprintf('%s/gab/',  cheminr);
index1 = 1;
chemin{2} = sprintf('%s/rgab/', cheminr);
index2 = 1;
chemin{3} = sprintf('%s/sgab', cheminr);
index3 = 1;
chemin{4} = sprintf('%s/vsgab', cheminr);
index4 = 1;

for i = 1:length(chemin)
 if ~exist(chemin{i}, 'dir')
   mkdir(chemin{i});
 end
end

gamma_K = 10^(-4);
f      = ((1:M2)-1)/M * Fs;
nfreqs = f / Fs;

fp = fopen(sprintf('%s/results_Fig2.txt',cheminr), 'w+b');

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
     
  T = w0T_range(idx1)/w0;
  L = T/Ts;

  tmp = sprintf('++ L=%0.02f, T=%0.02f ++\n', L, T);
  fprintf(1, tmp); fprintf(fp, tmp);


 %% 1) Compute the Gabor transform and display spectrogram
  %[tfrg, rtfrg] = tfrrgab(s, 1:length(s), M, round(L));

  [tfrg, rtfrg] = tfrrgab(s, M, L, gamma_K);   %% computed simultaneously with reassigned Gabor spectrogram
  
  %[tfrg] = tfrgab(s, M, L, gamma_K);  
  figure(1)
  imagesc(n, nfreqs, abs(tfrg(1:M2,:).^2).^alpha);
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Spectrogram SNR=%0.02f dB, L=%0.02f', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/gab_%d.eps', chemin{1}, index1), 'epsc');
  index1 = index1+1;
  
  %% 1b) reconstruction
  [ s_hat ] = rectfrgab(tfrg, L, M);
  
  tmp = sprintf('Gabor transform reconstruction, RQF=%0.02f, T=%0.02f\n', RQF(s,s_hat'), T);
  fprintf(1, tmp); fprintf(fp, tmp);
  
  
  %% 2) Reassigned Gabor Transform
  figure(2)
  imagesc(n, nfreqs, rtfrg(1:M2,:).^alpha);
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Reassigned spectrogram SNR=%0.02f dB, L=%d', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/rgab_%d.eps', chemin{2}, index2), 'epsc');
  index2 = index2+1;
  
  
  %% 3) Synchrosqueezed Gabor transform

  [~, stfrg] = tfrsgab(s, M, L, gamma_K);
  figure(3)
  imagesc(n, nfreqs, abs(stfrg(1:M2,:).^2).^alpha);
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Squared modulus of the Synch. STFT SNR=%0.02f dB, L=%d', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/sgab_%d.eps', chemin{3}, index3), 'epsc');
  index3 = index3+1;
  
  %% signal reconstruction
  [ s_hat1 ] = rectfrsgab(stfrg, L, M);
  
  %% 
  tmp = sprintf('Synchrosqueezed Gabor transform reconstruction, RQF=%0.02f, T=%0.02f\n', RQF(s,s_hat1'), T);
  fprintf(1, tmp); fprintf(fp, tmp);
  
  
  %% 4) Synchrosqueezed Gabor transform
   [~, stfrg] = tfrvsgab(s, M, L, gamma_K);
   
   figure(4)
   imagesc(n, nfreqs, abs(stfrg(1:M2,:).^2).^alpha);
   ylim([0 0.499])
   set(gca,'YDir','normal')
   xlabel('time samples', 'FontSize', 16)
   ylabel('normalized frequency', 'FontSize', 16)
   title(sprintf('Squared modulus of the Vertical Synch. STFT SNR=%0.02f dB, L=%d', rsb_target, L), 'FontSize', 14);
   colormap gray;
   cmap = colormap;
   cmap = flipud(cmap);
   colormap(cmap);
   saveas(gcf, sprintf('%s/sgab_%d.eps', chemin{4}, index4), 'epsc');
   index4 = index4+1;
   
   %% signal reconstruction
   [ s_hat4 ] = rectfrsgab(stfrg, L, M);
  
   tmp = sprintf('Vertical Synchrosqueezed Gabor transform reconstruction, RQF=%0.02f, T=%0.02f\n', RQF(s,s_hat4'), T);
   fprintf(1, tmp); fprintf(fp, tmp);
  
  
  for idx2 = 1:length(mu_range)
   mu = mu_range(idx2);
   
   %% 2b) LM-Reassigned Gabor Transform
   [~, rtfrg, lost] = tfrlmrgab(s, M, L, mu, gamma_K);
   
   %[tfrg, rtfrg] = tfrlmrgab(s, 1:length(s), M, mu, round(L*7+1));
   figure(2)
   imagesc(n, nfreqs, rtfrg(1:M2,:).^alpha);
   ylim([0 0.499])
   set(gca,'YDir','normal')
   xlabel('time samples', 'FontSize', 16)
   ylabel('normalized frequency', 'FontSize', 16)
   title(sprintf('LM-reassigned spectrogram SNR=%0.02f dB, L=%d \\mu=%0.02f', rsb_target, L, mu), 'FontSize', 14);
   colormap gray;
   cmap = colormap;
   cmap = flipud(cmap);
   colormap(cmap);
   saveas(gcf, sprintf('%s/rgab_%d.eps', chemin{2}, index2), 'epsc');
   index2 = index2+1;
    
   %% 3b) LM-Synchrosqueezed Gabor Transform
   [~, stfrg] = tfrlmsgab(s, M, L, mu, gamma_K);
   figure(3)
   imagesc(n, nfreqs, abs(stfrg(1:M2,:).^2).^alpha);
   ylim([0 0.499])
   set(gca,'YDir','normal')
   xlabel('time samples', 'FontSize', 16)
   ylabel('normalized frequency', 'FontSize', 16)
   title(sprintf('Squared modulus of the LM-Synch. STFT SNR=%0.02f dB, L=%d, \\mu=%0.02f', rsb_target, L, mu), 'FontSize', 14);
   colormap gray;
   cmap = colormap;
   cmap = flipud(cmap);
   colormap(cmap);
   saveas(gcf, sprintf('%s/sgab_%d.eps', chemin{3}, index3), 'epsc');
   index3 = index3+1;
   %% signal reconstruction
   [ s_hat2 ] = rectfrsgab(stfrg, L, M);
  
   tmp = sprintf('Synchrosqueezed Gabor transform reconstruction, RQF=%0.02f, T=%0.02f, mu=%.2f\n', RQF(s,s_hat2'), T, mu);
   fprintf(1, tmp); fprintf(fp, tmp);
  end
  
  
 end    %% real end
end  %% SNR

fclose all;
for i = 1:length(chemin)
 eps2pdf(chemin{i}, 1);
end