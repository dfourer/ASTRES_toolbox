%  ASTRES TOOLBOX test
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Mar. 2016
clear all
close all

chemin = './figs';
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

addpath('RLMR');  %% recursive STFT and vertical synchrosqueezed Gabor transform
%addpath('SSA');  %% SSA
addpath('ST');    %% ST

selected_signal = 1;  %%choose the signal to process: can be 1 or 2
gamma_K = 10^(-3);
recompute_all = false;
mat_file = 'Ligo_compare.mat';

if recompute_all || ~exist(mat_file, 'file')


if selected_signal == 1
    load fig1-observed-Livingston.txt
    t=fig1_observed_Livingston(:,1);
    s=fig1_observed_Livingston(:,2);
    Ts=t(2)-t(1);
    Fs=1.0/Ts;
    figs_folder = 'figs_signal1';
elseif selected_signal == 2

    load fig1-observed-Hanford.txt
    t=fig1_observed_Hanford(:,1);
    s=fig1_observed_Hanford(:,2);
    Ts=t(2)-t(1);
    Fs=1.0/Ts;
    figs_folder = 'figs_signal2';
else
    error('Unknown signal')
end

if ~exist(figs_folder, 'dir')
 mkdir(figs_folder);
end
NbPoints=length(s);


%% TF analysis
w0T = 7;
f0T = w0T / (2*pi);

w0  = 2*pi*50;
f0  = w0 / (2*pi);
T   = w0T / w0;
L = round(T / Ts);
%as_range = [w0/(2*pi*eps) w0/(pi/Ts)];    %% [40Hz -- Fs/2Hz]


M=3000; %5000     %% number of frequency bins
mi=floor(20*M/Fs); mf=ceil(520*M/Fs);   %% frequency range [20Hz - 520Hz]
as_range = [w0/(2*pi*5) w0/(pi/Ts)];
af_range = 1./as_range;

% center signal
%% s = s-mean(s);

q_method    = 2;
if_method   = 1;
q_threshold = gamma_K;
mm = m_axis(M);
nfreqs = mm/M;

m_range = mi:mf;
is_freq = 1;
n = 1:length(s);
for m = 1:3
 if m == 1
   %% Gabor
   fprintf(1,'Gabor transform - reassignment\n');
   tic
   [tfr{1}, rtfr{1}] = tfrrgab(s, M, L, gamma_K);
   %Nh = round(2 * L * sqrt(2*log(1/gamma_K)));
   %[tfr{1},rtfr{1}] = tfrrgab(s,1:length(s), M, Nh);
   toc
   
   fprintf(1,'synchrosqueezed Gabor transform\n');
   tic
   [~, stfr{1}] = tfrsgab(s, M, L, gamma_K);
   toc
   
   fprintf(1,'vertically synchrosqueezed Gabor transform\n');
   tic
   [~, vstfr{1}] = tfrvsgab(s, M, L, q_method, if_method, gamma_K, q_threshold);
   toc
 elseif  m == 2
   %% CWT
   fprintf(1,'Morlet WT - reassignment\n'); 
   tic
   [tfr{2},  as, rtfr{2}] = rMW( s, M, Ts, T, w0, af_range, gamma_K, is_freq);
   toc
   fprintf(1,'synchrosqueezed Morlet WT\n');
   tic
   [~,  ~, stfr{2}] = sMW( s, M, Ts, T, w0, af_range, gamma_K, is_freq);
   toc
   
   fprintf(1,'vertically synchrosqueezed Morlet WT\n');
   tic
   [~,  ~, vstfr{2}] = vsMW( s, M, Ts, T, w0, af_range, gamma_K, is_freq, 1);
   toc
 elseif  m == 3
   %% ST
   fprintf(1,'S-transform - reassignment\n');
   tic
   [tfr{3}, rtfr{3}] = tfrrst(s, M, w0T, gamma_K);
   toc
   
   tic
   fprintf(1,'S-transform - synchrosqueezed\n');
   [~, stfr{3}] = tfrsst(s, M, w0T, gamma_K); 
   toc

   tic
   fprintf(1,'S-transform - vertical synchrosqueezing\n');
   [~, vstfr{3}] = tfrvsst(s, M, w0T, gamma_K);
   toc 
 end
end

%as = a_axis(M,  as_range, 1);
af = Ts * w0 ./ (2*pi*as);
[~, i1] = min(abs(af-nfreqs(mi)));
[~, i2] = min(abs(af-nfreqs(mf)));
m_range2 = i1:i2;
save(mat_file);
else
    load(mat_file);
end

% %% TFR compression factor (modify display contrast)
alpha=0.5;
alpha2=0.2;
%% Gabor
figure(1)
imagesc(n, nfreqs(m_range), abs(tfr{1}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrgab(tfr{1}, L, M, mm);
q = RQF(s.', real(s_hat));
title(sprintf('spectrogram, T=%0.4f, RQF=%.2f dB', T, q), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1a.eps', chemin), 'epsc');

figure(2)
imagesc(n, nfreqs(m_range), abs(stfr{1}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrsgab(stfr{1}, L, M);
q1 = RQF(s.', real(s_hat));
title(sprintf('squared modulus of the synch. STFT, T=%0.4f, RQF=%.2f dB', T,q1), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1b.eps', chemin), 'epsc');

figure(3)
imagesc(n, nfreqs(m_range), abs(vstfr{1}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrsgab(vstfr{1}, L, M);
q2 = RQF(s.', real(s_hat));
title(sprintf('squared modulus of the vertical synch. STFT, T=%0.4f, RQF=%.2f dB', T,q2), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1c.eps', chemin), 'epsc');

figure(4)
imagesc(n, nfreqs(m_range), abs(rtfr{1}(m_range,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('reassigned spectrogram, T=%0.4f', T), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1d.eps', chemin), 'epsc');

%% MW
figure(5)
imagesc(n, af(m_range2), abs(tfr{2}(m_range2,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
s_hat = 2*recMW(tfr{2}, w0, T, as);
q3 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('Morlet scalogram, f0T=%.2f, RQF=%.2f dB', f0T, q3), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1e.eps', chemin), 'epsc');

figure(6)
imagesc(n, af(m_range2), abs(stfr{2}(m_range2,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
s_hat = recsMW(stfr{2}, w0, T);
q4 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('squared modulus of the synch. MW, f0T=%.2f, RQF=%.2f dB', f0T, q4), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1f.eps', chemin), 'epsc');

figure(7)
imagesc(n, af(m_range2), abs(vstfr{2}(m_range2,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
s_hat = recsMW(vstfr{2}, w0, T);
q5 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('squared modulus of the vert. synch. MW, f0T=%.2f, RQF=%.2f dB', f0T, q5), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1g.eps', chemin), 'epsc');

figure(8)
imagesc(n, af(m_range2), abs(rtfr{2}(m_range2,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('reassigned scalogram, f0T=%.2f', f0T), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1h.eps', chemin), 'epsc');

%% ST
figure(9)
imagesc(n, nfreqs(m_range), abs(tfr{3}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrst(tfr{3}, mm, M);
q6 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('Stockwellogram,f0T=%0.4f, RQF=%.2f dB', f0T, q6), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1i.eps', chemin), 'epsc');

figure(10)
imagesc(n, nfreqs(m_range), abs(stfr{3}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrsst(stfr{3}, w0T, mm);
q7 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('squared modulus of the synch. S-transform, f0T=%0.4f, RQF=%.2f dB', f0T,q7), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1j.eps', chemin), 'epsc');

figure(11)
imagesc(n, nfreqs(m_range), abs(vstfr{3}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
[ s_hat ] = rectfrsst(vstfr{3}, w0T, mm);
q8 = RQF(s(:).', real(s_hat(:).'));
title(sprintf('squared modulus of the vertical synch. S-transform, f0T=%0.4f, RQF=%.2f dB', f0T,q8), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1k.eps', chemin), 'epsc');

figure(12)
imagesc(n, nfreqs(m_range), abs(rtfr{3}(m_range,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('reassigned Stockwellogram, f0T=%0.4f', f0T), 'FontSize', 10);
colormap gray; cmap = flipud(colormap);colormap(cmap);
saveas(gcf, sprintf('%s/Fig1l.eps', chemin), 'epsc');






figure(999)
subplot(331)
imagesc(n, nfreqs(m_range), abs(tfr{1}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('Gabor spectrogram, T=%0.4f', T), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(332)
imagesc(n, af(m_range2), abs(tfr{2}(m_range2,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('Morlet scalogram, f0T=%.2f', f0T), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(333)
imagesc(n, nfreqs(m_range), abs(tfr{3}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('Stockwellogram, f0T=%.2f', f0T), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(334)
imagesc(n, nfreqs(m_range), abs(stfr{1}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('squ. mod. synch. STFT, T=%0.4f, RQF=%.2f dB', T, q1), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(335)
imagesc(n, af(m_range2), abs(stfr{2}(m_range2,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('squ. mod. synch. MW, f0T=%.2f, RQF=%.2f dB', f0T, q4), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(336)
imagesc(n, nfreqs(m_range), abs(stfr{3}(m_range,:).^2).^alpha);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
title(sprintf('squ. mod. synch. ST, f0T=%.2f, RQF=%.2f dB', f0T, q7), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(337)
imagesc(n, nfreqs(m_range), abs(vstfr{1}(m_range,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
% [ s_hat ] = rectfrsgab(stfr{1}, L, M);
% rqf_stfr  = RQF(s.', real(s_hat));
title(sprintf('squ. mod. vertical synch. STFT T=%0.4f, RQF=%.2f dB', T, q2), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(338)
imagesc(n, af(m_range2), abs(vstfr{2}(m_range2,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
% [ s_hat ] = recsMW(stfr{2}, w0, T); %, m_range
% rqf_sWT  = RQF(s, real(s_hat));
title(sprintf('squ. mod. vertical synch. MW, f0T=%.2f, RQF=%.2f dB', f0T, q5), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

subplot(339)
imagesc(n, nfreqs(m_range), abs(vstfr{3}(m_range,:).^2).^alpha2);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 10)
ylabel('normalized frequency', 'FontSize', 10)
% [s_hat] = rectfrsst(stfr{3}, w0T); %, mm
% rqf_SST  = RQF(s.', real(s_hat));
title(sprintf('squ. mod. vertical synch. ST, f0T=%.2f, RQF=%.2f dB', f0T, q8), 'FontSize', 8);
colormap gray; cmap = flipud(colormap);colormap(cmap);

saveas(gcf, sprintf('%s/Fig_all.eps', chemin), 'epsc');




eps2pdf(chemin);