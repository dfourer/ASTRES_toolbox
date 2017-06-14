%%   Script used to generate all figures
%
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 24-07-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassigning and Synchrosqueezing of the Stockwell Transform. IEEE. Tran. on Signal Processing. 2015]
clear all
close all

signal = 1; %5;  %% change value in [1,2,3,4]

mu_range  = [0.01 0.05 0.2 0.5 1 2 5];      %mu_range = logspace(log10(mu_range(1)), log10(mu_range(2)), NbImg2);
SNR_range = [45 25 10];%[45 35 25 15 5];
%w0T_range = [2*pi*2 2*pi*5 2*pi*10];        %L_range = linspace(2*pi*1, 2*pi*15, NbImg); % (2*pi*10)          %% range for omega_0 T
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
s0 = real(s-mean(s)); %% 0-mean signal

%%% /!\ Use a linear scale to obtain a correct mapping using imagesc
w0 = 2*pi*100;
%as_range = [w0/(2*pi*eps) w0/(pi/Ts)];    %% [40Hz -- Fs/2Hz]

as_range = [w0/(2*pi*5) w0/(pi/Ts)];
%as_range = [1.4 0.2];
af_range = 1./as_range;   %% used to compute scales from normalized frequencies

%%% prepare folders to store images
chemin0 = 'figs_TSP';

cheminr = sprintf('%s/%d', chemin0, signal);
chemin{1} = sprintf('%s/mwt/',  cheminr);
index1 = 1;
chemin{2} = sprintf('%s/rmwt/', cheminr);
index2 = 1;
chemin{3} = sprintf('%s/smwt', cheminr);
index3 = 1;


for i = 1:length(chemin)
 if ~exist(chemin{i}, 'dir')
   mkdir(chemin{i});
 end
end

gamma_K = 10^(-4);
f      = ((1:M2)-1)/M * Fs;
nfreqs = f / Fs;

fp = fopen(sprintf('%s/results_Fig3.txt',cheminr), 'w+b');

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
  
  T = w0T / w0;
  L = T / Ts;
  
  tmp = sprintf('++ L=%0.02f, T=%0.02f ++\n', L, T);
  fprintf(1, tmp); fprintf(fp, tmp);
  
  
  %% 1) Morlet wavelet transform
  %[Wx,as,dWx] = cwt_fw(s, 'morlet', nv, Ts,  struct('mu', L));
  [ Wx, as ] = MW( s, M, Ts, T, w0, af_range, gamma_K, 1);
  af = Ts * w0 ./ (2*pi*as);
  
  figure(1)
  imagesc(n, af, (abs(Wx).^2).^alpha);
 %imagesc(n, af, (abs(stfr).^2).^alpha);
  %set(gca,'Ytick',0:0.05:0.5)
  %set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Morlet scalogram SNR=%0.02f dB, L=%0.2f', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  
% [ Wx, as ] = MW( s, M, Ts, T, w0, as_range, gamma_K, 2);
%   figure(1)
%   imagesc(n, as, (abs(Wx).^2).^alpha);
%   %set(gca,'YDir','normal')
%   xlabel('time samples', 'FontSize', 16)
%   ylabel('scale', 'FontSize', 16)
%   title(sprintf('Morlet scalogram SNR=%0.02f dB, L=%0.2f', rsb_target, L), 'FontSize', 14);
%   colormap gray;
%   cmap = colormap;
%   cmap = flipud(cmap);
%   colormap(cmap);
  saveas(gcf, sprintf('%s/mwt_%d.eps', chemin{1}, index1), 'epsc');
  index1 = index1+1;

  %% 3.1) Reconstruction 
  s_hat = real(recMW(Wx, w0, T, as));
  
%   plot(s)
%   hold on
%   plot(s_hat, 'r-.')
   res = RQF(s, s_hat');
  
  tmp = sprintf('Morlet Wavelet Transform reconstruction (Truncated Morlet formula), L=%0.2f, T=%0.2f, f0=%0.2f, RQF=%0.2f \n ', L, T, w0/(2*pi), res);
  fprintf(1, tmp); fprintf(fp, tmp);
  
  
  
  %% 2) Reassigned Morlet wavelet transform
  [Wx,  as, rtfr, af, lost] = rMW( s, M, Ts, T, w0, af_range, gamma_K, 1);
  %[ Wx, as, rtfr, af, lost ] = rMW( s, M, Ts, T, w0, as_range, gamma_K );
  figure(2)
  imagesc(n, af*Ts, rtfr.^alpha);  %af*Ts
%   set(gca,'Ytick',0:0.05:0.5)
%   set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  %ylabel('scale', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Reassigned Morlet scalogram SNR=%0.02f dB, L=%0.2f', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/rmwt_%d.eps', chemin{2}, index2), 'epsc');
  index2 = index2+1;
  
  
  %% 3) synchrosqueezed Morlet wavelet transform
  [ Wx, as, stfr, af, lost ] = sMW( s, M, Ts, T, w0, af_range, gamma_K,1);
  figure(3)
  imagesc(n, af*Ts, (abs(stfr).^2).^alpha);  %af*Ts
%   set(gca,'Ytick',0:0.05:0.5)
%   set(gca,'Xtick',0:50:N)
  ylim([0 0.499])
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  %ylabel('scale', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Squared modulus of the synch. Morlet WT SNR=%0.02f dB, L=%0.2f', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  saveas(gcf, sprintf('%s/smwt_%d.eps', chemin{3}, index3), 'epsc');
  index3 = index3+1;
  
  %% 3.1) Reconstruction 
  s_hat = real(recsMW(stfr, w0, T));
  
%   plot(s)
%   hold on
%   plot(s_hat, 'r-.')
   res = RQF(s, s_hat);
  
  tmp = sprintf('Synchrosqueezed Morlet Wavelet Transform, L=%0.2f, T=%0.2f, f0=%0.2f, RQF=%0.2f \n ', L, T, w0/(2*pi), res);
  fprintf(1, tmp); fprintf(fp, tmp);
  
  %% 4) vertical Synchrosqueezed Morlet wavelet transform
  [ Wx, as, stfr, af, lost ] = vsMW( s, M, Ts, T, w0, af_range, gamma_K,1);
  figure(4)
  imagesc(n, af*Ts, (abs(stfr).^2).^alpha);  %af*Ts
  ylim([0 0.499]);
  set(gca,'YDir','normal')
  xlabel('time samples', 'FontSize', 16)
  %ylabel('scale', 'FontSize', 16)
  ylabel('normalized frequency', 'FontSize', 16)
  title(sprintf('Squared modulus of the synch. Morlet WT SNR=%0.02f dB, L=%0.2f', rsb_target, L), 'FontSize', 14);
  colormap gray;
  cmap = colormap;
  cmap = flipud(cmap);
  colormap(cmap);
  %saveas(gcf, sprintf('%s/smwt_%d.eps', chemin{3}, index3), 'epsc');
  %index3 = index3+1;
  
  %% 3.1) Reconstruction 
  s_hat = real(recsMW(stfr, w0, T));
%   plot(s)
%   hold on
%   plot(s_hat, 'r-.')
   res = RQF(s, s_hat);
  
  tmp = sprintf('Synchrosqueezed Morlet Wavelet Transform, L=%0.2f, T=%0.2f, f0=%0.2f, RQF=%0.2f \n ', L, T, w0/(2*pi), res);
  fprintf(1, tmp);
  %fprintf(fp, tmp);
  
  
%   w = cwt_freq(Wx, Ts);
%   [Tx,fs] = synsq_cwt_squeeze(Wx, w, t, nv);
%   figure(3)
%   imagesc(n, fs, (abs(Tx).^2).^alpha);
%   set(gca,'YDir','normal')
%   xlabel('time samples', 'FontSize', 16)
%   ylabel('frequency [Hz]', 'FontSize', 16)
%   title(sprintf('square modulus of the synchrosqueezed CWT SNR=%0.02f dB, L=%0.02f', rsb_target, L), 'FontSize', 14);
%   colormap gray;
%   cmap = colormap;
%   cmap = flipud(cmap);
%   colormap(cmap);
%   saveas(gcf, sprintf('%s/mwt_%d.eps', chemin{3}, index3), 'epsc');
%   index3 = index3+1;
  
 end    %% real end
end  %% SNR

fclose all;
for i = 1:length(chemin)
 eps2pdf(chemin{i}, 1);
end