%% test non recursive vertical synchrosqueezing

clear all
close all

if_method = 2;

signal=1;
%% used for zero-padding for Reconstruction Quality computation
N_0pad = 100;
z      = zeros(N_0pad, 1);


%% 1: Load test signal
M = 600;
mi = 1;mf = round(M/2);

k_range   = [5 7];
L_range   = [5 7 10];
SNR_range = [45 25 10];
%n0_range0 = [8 18 27];
%mu_range = [0.3 0.8 1.3 1.8 2.3 2.8];
M_range = [100 200 600 1000 2400];

alpha  = 0.3;           %% used for display

load_signal;


s = sigmerge(s,randn(size(s)), SNR_range(1));


%% remove the first impulse
% s(1:60) = 0;
% s(30) = 20;

t = 1:N;
n_freq = m_axis(M)/M;
Mh = round(M/2);
m_range = 1:Mh;

%% Vertical synchrosqueezing using CRE1
L = 8;  %20
q_method = 1;

gamma_K = 1e-4;
[tfr, stfr, lost, q_hat_map] = my_tfrvsgab(s, M, L, q_method, if_method, gamma_K);
[ s_hat ] = my_rectfrsgab(stfr, L, M);

figure(1)
imagesc(t, n_freq(m_range), abs(stfr(m_range,:)).^alpha);
set(gca,'YDir','normal')
colormap gray; cmap = colormap;
colormap(flipud(cmap)); 
xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16);
title(sprintf('Squared modulus of the Vertical synchrosqueezed STFT q_method=%d, if_method=%d', q_method, if_method));

%% Vertical synchrosqueezing using CRE2t
q_method = 2;
[tfr2, stfr2, lost2, q_hat_map2] = my_tfrvsgab(s, M, L, q_method, if_method, gamma_K);
[ s_hat2 ] = my_rectfrsgab(stfr2, L, M);

figure(2)
imagesc(t, n_freq(m_range), abs(stfr2(m_range,:)).^alpha);
set(gca,'YDir','normal')
colormap gray; cmap = colormap;
colormap(flipud(cmap)); 
xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16);
title(sprintf('Squared modulus of the Vertical synchrosqueezed STFT q_method=%d, if_method=%d', q_method, if_method));

% figure(3)
% imagesc(t, nfreqs, 10 * log10(abs(q_hat_map)));
% title(sprintf('Area where q(t,\\omega) is computed (white area) is computed k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
% 
% figure(4)
% imagesc(t, nfreqs, 10 * log10(abs(q_hat_map2)));
% title(sprintf('Area where q(t,\\omega) is computed (white area) is computed k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);

%% Vertical synchrosqueezing using CRE2w
q_method = 3;
[tfr3, stfr3, lost3, q_hat_map3] = my_tfrvsgab(s, M, L, q_method, if_method, gamma_K);
[ s_hat3 ] = my_rectfrsgab(stfr3, L, M);

figure(3)
imagesc(t, n_freq(m_range), abs(stfr3(m_range,:)).^alpha);
set(gca,'YDir','normal')
colormap gray; cmap = colormap;
colormap(flipud(cmap)); 
xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16);
title(sprintf('Squared modulus of the Vertical synchrosqueezed STFT q_method=%d, if_method=%d', q_method, if_method));

%% Vertical synchrosqueezing using CRE2r
q_method = 4;
[tfr4, stfr4, lost4, q_hat_map4] = my_tfrvsgab(s, M, L, q_method, if_method, gamma_K);
[ s_hat4 ] = my_rectfrsgab(stfr4, L, M);

figure(4)
imagesc(t, n_freq(m_range), abs(stfr4(m_range,:)).^alpha);
set(gca,'YDir','normal')
colormap gray; cmap = colormap;
colormap(flipud(cmap)); 
xlabel('time samples', 'FontSize', 16);  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16);
title(sprintf('Squared modulus of the Vertical synchrosqueezed STFT q_method=%d, if_method=%d', q_method, if_method));


