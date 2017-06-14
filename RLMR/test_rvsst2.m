%% comparison between chirp-rate estimator (using recursive VSST)

clear all
close all


chemin = 'figs';    %% reass
if ~exist(chemin, 'dir')
 mkdir(chemin);
end


%% used for zero-padding for Reconstruction Quality computation
N_0pad = 100;
z      = zeros(N_0pad, 1);


M   = 600;
mi  = 1; mf = round(M/2);
k   = 5;
L   = 7;
n0 = (k-1) * L;
contrst  = 0.35;           %% used for display


%%%    Gaussian modulated Amplitude chirp + sinusoid
N = 500;
alpha = 0.35*2*pi/N;
t = (0:N-1);
A   = 1;
t0 = 250;
T_x = 200; %inf; %200;   %% set T_x to inf for constant amplitude
Ax    = A * exp(-(t-t0).^2 / (2*T_x^2));
phi_x =  alpha * t.^2/2;  %2*pi * 0 *t +

s = real(Ax .* exp(1j * phi_x)) + real(Ax .* exp(1j .* 2*pi * 0.4 * t));
s = s(:);

n_range0 = N_0pad+(1:N);

s_ref = s;
Mr = mf-mi+1;

rsb_target = 50;
s = sigmerge(s_ref,randn(size(s_ref)),rsb_target);
s = [z;s];                         %% 0 padding used for reconstruction
rsb = SNR(s_ref, s(n_range0));


%% classical synchrosqueezing
[tfr, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0);

figure(1) 
imagesc(t, nfreqs, (abs(tfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Recurs. spectrogram k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/spectrogram.eps', chemin), 'epsc');


figure(2) 
imagesc(t, nfreqs, (abs(stfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. synchrosqueezed STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rsst.eps', chemin), 'epsc');

%% reconstruction Sychrosqueezing
[ s_hat] = sstft_rec(stfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Classical synchrosqueezing reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);


%% vertical synchrosqueezing  CRE1

%% CRE1 with biased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 1, 1);

figure(31)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE1 - IF biased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE1_IFbiased.eps', chemin), 'epsc');


% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE1 - IF biased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);


%% CRE1 with unbiased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 1, 2);
figure(32)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE1 - IF unbiased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE1_IFunbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE1 - IF unbiased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);


%% CRE2_t with biased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map_cre2t] = recursive_vsstft(s, k, L, mi, mf, M, n0, 2, 1);
figure(41)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE2_t - IF biased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE2t_IFbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE2_t - IF biased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);

%% CRE2_t with unbiased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 2, 2);
figure(42)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE2_t - IF unbiased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE2t_IFunbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE2_t - IF unbiased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);


%% CRE2_w with biased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map_cre2w] = recursive_vsstft(s, k, L, mi, mf, M, n0, 3, 1);
figure(51)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE2_w - IF biased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE2w_IFbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE2_w - IF biased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);

%% CRE2_w with unbiased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 3, 2);
figure(52)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE2_w - IF unbiased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE2w_IFunbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE2_w - IF unbiased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);


%% CRE2_robust with biased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map_crerobust] = recursive_vsstft(s, k, L, mi, mf, M, n0, 4, 1);
figure(61)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE robust - IF biased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE_robust_IFbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE_robust - IF biased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);

%% CRE2_robust with unbiased IF estimator
[~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 4, 2);
figure(62)
imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^contrst);
set(gca,'YDir','normal')
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Squared modulus of the Recurs. vert. synchrosqueezed (CRE robust - IF unbiased) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
saveas(gcf, sprintf('%s/rvsst_CRE_robust_IFunbiased.eps', chemin), 'epsc');

% reconstruction
[ s_hat] = sstft_rec(vstfr, k, L, M, n0);
n_range_rec = n_range0-n0;
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
tmp_res = sprintf('Vertical synchrosqueezing (CRE_robust - IF unbiased) reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
fprintf(1,  tmp_res);





eps2pdf(chemin);








