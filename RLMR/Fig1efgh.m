%
%  Generate Fig 1 (e)-(h) of Ref
%  Performance curves of Chirp Rate (CR) and Instantaneous Frequency (IF)
%  estimators
%
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 27-06-2016
% Ref: [D. Fourer, F. Auger, K.Czarnecki, S. Meignen and P. Flandrin, Chirp rate and instantaneous frequency estimation. Proc. IEEE ICASSP 2017]

clear all
close all

warning ('off','all');
grid_enable = true;

chemin = 'figs';
if ~exist(chemin, 'dir')
 mkdir(chemin);
end

use_recursive = true;  % use the recursive implementation or not

signal=12;
%% used for zero-padding for Reconstruction Quality computation
%N_0pad = 100;
%z      = zeros(N_0pad, 1);

q_methods  = {'CRE1' 'CRE2_t' 'CRE2_\omega' 'CRE2_r'};
if_methods = {'IF(1)' 'IF(2)'};

nb_qmethods  = 4;
nb_ifmethods = 2;
%% 1: Load test signal
M = 500;
Mh = round(M/2);

if use_recursive
 k = 7; %7
 L = 6; %6;
 n0 = (k-1) * L;
 n0r = 75;
 mi = 0;
 mf = Mh;
 q_threshold = 1e-4; %eps; %1e-2;
else
 L = 20;   
 q_threshold = 1e-4;
end

gamma_K = 1e-4;

SNR_range = [45 25 10]; %[45 25 10];
alpha  = 0.4; %0.3           %% used for display

%mi = 1;mf = round(M/2);
%k_range   = [5 7];
%L_range   = [5 7 10];
%n0_range0 = [8 18 27];
%mu_range = [0.3 0.8 1.3 1.8 2.3 2.8];
%M_range = [100 200 600 1000 2400];
load_signal;

t = 1:N;
n_freq = m_axis(M)/M;
Mh = round(M/2);
m_range = 1:Mh;

index = ones(1,4);
for i_snr = 1:length(SNR_range)

  snr = SNR_range(i_snr);
  x = sigmerge(s,randn(size(s)), snr);

  
  if use_recursive
   [tfr, stfr, nfreqs, below, over] = recursive_sstft(x, k, L, mi, mf, M, n0);
   titre1 = sprintf('$|F_x^h(t,\\omega)|^2$, $k$=%d, $L_h$=%d, SNR=%d dB', k, L,round(snr));
   titre2 = sprintf('$|SF_x^h(t,\\omega)|^2$, $k$=%d, $L_h$=%d, SNR=%d dB', k, L, round(snr));
  else
   [tfr, stfr, lost] = my_tfrsgab(x, M, L, gamma_K);
   titre1 = sprintf('$|F_x^h(t,\\omega)|^2$, $L_h$=%d, SNR=%d dB', L, round(snr));
   titre2 = sprintf('$|SF_x^h(t,\\omega)|^2$, $L_h$=%d, SNR=%d dB', L, round(snr));
  end

  %% 1 classical spectrogram
  figure(1)
  %% reconstruction quality
  if use_recursive
   [ s_hat ] = stft_rec(tfr, k, L, mi, mf, M, n0);
   n_range_rec = n0r:(length(s_hat));
   s_rec = 2*real(s_hat(n_range_rec));    %% assume the signal is real
  else
    [ s_hat ] = my_rectfrsgab(tfr, L, M);
    s_rec = real(s_hat(n_range_rec));
    n_range_rec = 1:length(s_rec);
  end
  rq = RQF(x(n_range_rec), s_rec);
  imagesc(t, n_freq(m_range), (abs(tfr(m_range,:)).^2).^alpha);
  if grid_enable,   grid on; end
  set(gca,'YDir','normal')
  colormap gray; cmap = colormap;
  colormap(flipud(cmap));
  title(sprintf('%s, RQF=%.2f dB',titre1, rq), 'FontSize', 15, 'Interpreter','Latex');
  xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
  ylabel('normalized frequency', 'FontSize', 16);
  saveas(gcf, sprintf('%s/fig1e_spectrogram_%d.eps', chemin, index(1)), 'epsc');
  index(1) = index(1)+1;
  

  %% 1 sq. mod. classical synchrosqueezed STFT
  figure(2)
  %% reconstruction quality
  if use_recursive
   [ s_hat] = sstft_rec(stfr, k, L, M, n0);
   n_range_rec = n0r:(length(s_hat));
   s_rec = 2*real(s_hat(n_range_rec));    %% assume the signal is real
  else
    [ s_hat ] = my_rectfrsgab(stfr, L, M);
    s_rec = real(s_hat(n_range_rec));
    n_range_rec = 1:length(s_rec);
  end
   rq = RQF(x(n_range_rec), s_rec);
   
  imagesc(t, n_freq(m_range), (abs(stfr(m_range,:)).^2).^alpha);
  if grid_enable,   grid on; end
  set(gca,'YDir','normal')
  colormap gray; cmap = colormap;
  colormap(flipud(cmap));
  title(sprintf('%s, RQF=%.2f dB', titre2,rq), 'FontSize', 15, 'Interpreter','Latex');
  xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
  ylabel('normalized frequency', 'FontSize', 16); 
  saveas(gcf, sprintf('%s/fig1f_sstft_%d.eps', chemin, index(2)), 'epsc');
  index(2) = index(2)+1;
  
  %% Vertical synchrosqueezing using CRE1
  ifq_methods = {'$\hat{\dot{\phi}}_x^{K1}$' '$\hat{\dot{\phi}}_x^{t2}$' '$\hat{\dot{\phi}}_x^{\omega 2}$' '$\hat{\dot{\phi}}_x^{r2}$' '$\hat{\omega}$'};
  
  q_method = 1;   if_method = 1;
   if use_recursive
    [~, vstfr] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_method, if_method, q_threshold);
    titre3 = sprintf('$|VSF_x^h(t,\\omega)|^2$ using %s, $k$=%d, $L_h$=%d, SNR=%d dB', ifq_methods{1}, k, L, round(snr));
   else
    [~, vstfr] = my_tfrvsgab(x, M, L, q_method, if_method, gamma_K, q_threshold);
    titre3 = sprintf('$|VSF_x^h(t,\\omega)|^2$ using %s, $L_h$=%d, SNR=%d dB', ifq_methods{1}, L, round(snr));
   end
   
   %[ s_hat ] = my_rectfrsgab(stfr, L, M);
   figure(3)
   if use_recursive
    [ s_hat] = sstft_rec(vstfr, k, L, M, n0);
    n_range_rec = n0r:(length(s_hat));
    s_rec = 2*real(s_hat(n_range_rec));    %% assume the signal is real
   else
    [ s_hat ] = my_rectfrsgab(vstfr, L, M);
    s_rec = real(s_hat(n_range_rec));
    n_range_rec = 1:length(s_rec);
    end
   rq = RQF(x(n_range_rec), s_rec);
   
   imagesc(t, n_freq(m_range), (abs(vstfr(m_range,:)).^2).^alpha);
   if grid_enable,   grid on; end
   set(gca,'YDir','normal')
   colormap gray; cmap = colormap;
   colormap(flipud(cmap));
   title(sprintf('%s, RQF=%.2f dB', titre3, rq), 'FontSize', 15, 'Interpreter','Latex');
   xlabel('time samples', 'FontSize', 16);
   ylabel('normalized frequency', 'FontSize', 16);
   saveas(gcf, sprintf('%s/fig1g_vsstft_%d.eps', chemin, index(3)), 'epsc');
   index(3) = index(3)+1;
 
  
  %% Vertical synchrosqueezing using CRE2
  for q_method = 2:length(q_methods)
    for if_method = 2
        
     if use_recursive
      [~, vstfr2] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_method, if_method, q_threshold);
      titre4 = sprintf('$|VSF_x^h(t,\\omega)|^2$ using %s, $k$=%d, $L_h$=%d, SNR=%d dB', ifq_methods{q_method}, k, L,round(snr));
     else
      [~, vstfr2] = my_tfrvsgab(x, M, L, q_method, if_method, gamma_K, q_threshold);
      titre4 = sprintf('$|VSF_x^h(t,\\omega)|^2$ using %s, $L_h$=%d, SNR=%d dB', ifq_methods{q_method}, L, round(snr));
     end
     
     figure(4)
     if use_recursive
      [ s_hat] = sstft_rec(vstfr2, k, L, M, n0);
      n_range_rec = n0r:(length(s_hat));
      s_rec = 2*real(s_hat(n_range_rec));    %% assume the signal is real
     else 
      [ s_hat ] = my_rectfrsgab(vstfr2, L, M);
      s_rec = real(s_hat(n_range_rec));
      n_range_rec = 1:length(s_rec);
     end 
     rq = RQF(x(n_range_rec), s_rec);
     imagesc(t, n_freq(m_range), (abs(vstfr2(m_range,:)).^2).^alpha);
     if grid_enable,   grid on; end
     set(gca,'YDir','normal')
     colormap gray; cmap = colormap;
     colormap(flipud(cmap)); 
     title(sprintf('%s, RQF=%.2f dB', titre4, rq), 'FontSize', 15, 'Interpreter','Latex');
     %title(sprintf('Squared modulus of the vertical synchrosqueezed STFT (%s - %s), SNR=%d, L_h=%d', q_methods{q_method}, if_methods{if_method}, round(snr), L), 'FontSize', 12);
     xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
     ylabel('normalized frequency', 'FontSize', 16);
     saveas(gcf, sprintf('%s/fig1h_vsstft2_%d.eps', chemin, index(4)), 'epsc');
     index(4) = index(4)+1;
    end
  end
end
eps2pdf(chemin);
warning ('on','all');
