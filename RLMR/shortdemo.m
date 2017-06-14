% Recursive versions of the Levenberg-Marquardt reassigned spectrogram 
% and of the synchrosqueezed STFT short demo
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Feb. 2016
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]
clear all
close all

alpha=0.5;             %% TFR compression factor (modify display contrast)
N_0pad = 50;           %% number of zeros to add before signal (used for reconstruction)
AddedZeros = zeros(N_0pad, 1);

%% input signal (3 modes with zero-padding)
S = zeros(3, N_0pad+3*60+N_0pad);
S(1,N_0pad+(1:80))     = real(fmconst(80,0.1));
S(2,N_0pad+60+(1:80))  = real(fmconst(80,0.2));
S(3,N_0pad+120+(1:80)) = real(fmconst(80,0.3));

s = sum(S).';
N = length(s)-N_0pad;         %% signal length without 0-padding
t = 1:length(s);

figure(1);
plot(s);
xlabel('time samples', 'FontSize', 16)
ylabel('amplitude', 'FontSize', 16)
title('Analyzed signal')



%% TF analysis
fprintf(1,'Press a key to start the time-frequency analysis\n');
pause

k=3;                   %% recursive filter order
L=10;                  %% analysis window length
mi=1; mf=250; M=500;   %% frequency range

%% Recursive spectrogram
fprintf(1,'Computing the recursive respectrogram and the recursive reassigned respectrogram...\n');
[tfr, nfreqs] = recursive_stft(s, k, L, mi, mf, M);
[~, rtfr, nfreqs, before, after, below, over] = recursive_rsp(s, k, L, mi, mf, M);

figure(2);
imagesc(t, nfreqs, (abs(tfr).^2).^alpha);
set(gca,'YDir','normal')
colormap gray;
colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Recursive spectrogram k=%d, L=%2.2f ', k, L), 'FontSize', 14);

%% Recursive reassigned spectrogram
figure(3);
imagesc(t, nfreqs, rtfr.^alpha);
set(gca,'YDir','normal')
colormap gray;
colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Recursive reassigned spectrogram k=%d, L=%2.2f ', k, L), 'FontSize', 14);


%% Recursive Levenberg-Marquardt reassigned spectrogram
fprintf(1,'Press a key to compute the recursive Levenberg-Marquardt reassigned spectrogram\n');
pause
mu_range = 0.5:0.5:8;
for im = 1:length(mu_range)
 mu = mu_range(im);
 [~, rtfr, nfreqs, before, after, below, over] = recursive_lmrsp(s, k, L, mi, mf, M, mu);
 figure(31);
 imagesc(t, nfreqs, rtfr.^alpha);
 set(gca,'YDir','normal')
 colormap gray;
 colormap(flipud(colormap)); 
 xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
 ylabel('normalized frequency', 'FontSize', 16)
 title(sprintf('Recursive LM-reassigned spectrogram k=%d, L=%2.2f, mu=%.2f ', k, L, mu), 'FontSize', 14);
 F(im) = getframe;

end
title(sprintf('Recursive LM-reassigned spectrogram k=%d, L=%.2f,mu in [%.2f,%.2f]', k, L, mu_range(1), mu_range(end)), 'FontSize', 14);
movie(F, 3); %%replay 3 times


fprintf(1,'Press a key to compute the recursive synchrosqueezed STFT\n');
pause
%% Recursive synchrosqueezed STFT
n0 = round((k-1) * L);

if n0 > N_0pad
 warning(0, 'N_0pad should be chosen greater than n0=%d', n0);
end

[~, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0);
figure(4);
imagesc(t, nfreqs, (abs(stfr).^2).^alpha);
set(gca,'YDir','normal')
colormap gray;
colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Recursive synchrosqueezed STFT k=%d, L=%2.2f ', k, L), 'FontSize', 14);

%% Signal reconstruction from recursive STFT
n_range     = N_0pad+(1:N);
n_range_rec = n_range-n0;
s_hat = stft_rec(tfr, k, L, mi, mf, M, n0);
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
rqf_stft = RQF(s(n_range_rec), s_rec);
fprintf(1,'Recursive STFT reconstruction : RQF=%0.02f dB\n', rqf_stft);

figure(11);
plot(t(n_range_rec), s(n_range_rec), 'g-.');
hold on
plot(t(n_range_rec), s_rec, 'k');
legend('original signal', 'reconstructed signal')
title('recursive STFT reconstruction')
title(sprintf('recursive synchrosqueezed STFT reconstruction, RQF=%.2f dB', rqf_stft))


%% Signal reconstruction from recursive synchrosqueezed STFT
s_hat = sstft_rec(stfr, k, L, M, n0);
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
rqf_sstft =  RQF(s(n_range_rec), s_rec);
fprintf(1,'Recursive synchrosqueezed STFT reconstruction : RQF=%0.02f dB\n', rqf_sstft);

figure(12);
plot(t(n_range_rec), s(n_range_rec), 'g-.');
hold on
plot(t(n_range_rec), s_rec, 'k');
legend('original signal', 'reconstructed signal')
title(sprintf('recursive synchrosqueezed STFT reconstruction, RQF=%.2f dB', rqf_sstft))


%% Mode separation from recursive synchrosqueezed STFT
mask_sine1 = zeros(size(stfr));
m_ridge = round(0.2*M-1); delta_m = 8;
mask_sine1((m_ridge-delta_m):(m_ridge+delta_m), 130:210) = 1;
    
filtered_stfr = stfr .* mask_sine1;
figure(41);
imagesc(abs(filtered_stfr));
set(gca,'YDir','normal')
colormap gray;
colormap(flipud(colormap)); 
xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
ylabel('normalized frequency', 'FontSize', 16)
title(sprintf('Filtered recursive synchrosqueezed STFT k=%d, L=%2.2f ', k, L), 'FontSize', 14);

s_hat = sstft_rec(filtered_stfr, k, L, M, n0);
s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
 
figure(13)
plot(s_rec)
hold on
plot(S(2,n_range_rec), 'g-.')
legend('reconstructed signal', 'original signal')
rqf_mode2 =  RQF(S(2,n_range_rec).', s_rec);
title(sprintf('reconstructed 2nd mode, RQF=%.2f', rqf_mode2))


