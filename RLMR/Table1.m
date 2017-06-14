% Example script for recursive STFT and recursive synchrosqueezed STFT
% reconstruction
% 
%
% select the test signal setting the variable signal in range 1-6
% signal=1; load_signal;
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

clear all
close all

%% used for zero-padding for Reconstruction Quality computation
N_0pad = 100;
z      = zeros(N_0pad, 1);

%% 1: Load test signal
M = 300;
mi = 1;mf = round(M/2);

plot_snr_curve = true;        %% to plot SNR curves

k_range   = 5; %[5 9];        %5:-1:3;   %% 3
L_range   = 7;            %8:-1:6;   %% 3
SNR_range = 45; %[45 35 25 15 5]; %[45 20 5];      %60:-20:0; %% 3

n0_range0 = [8 18 26];
mu_range = [0.3 0.8 1.3 1.8 2.3]; % 2.8]; %[0.45 0.6 0.8 1.25 2.4]; %10
M_range = [100 200 600 1000 2400];

signal=1; load_signal;
figure(6); 
subplot(211); plot(real(s));title('real(s)')
subplot(212); plot(imag(s));title('imag(s)')

n_range0 = N_0pad+(1:N);

s_ref = s;
Mr = mf-mi+1;

chemin  = sprintf('figs/%d/', signal);
chemin3 = sprintf('%s/table', chemin);

if ~exist(chemin3, 'dir')
 mkdir(chemin3);
end

fp = fopen(sprintf('%s/results.txt', chemin3), 'wt');

%% 2: Recursive reassignment/synchrosqueezing parameters
% M: number of frequency bins, mi,mf: frequencies bounds, k: order,
% L:window length

%Ts = 1.0 / Fs;
t = (1:length(s))-1;


%% 3.1) Compute the Recursive spectrogram
index  = 1;     %% used for stft
index2 = 1;     %% used for reassignment LM reassigned spectr
index3 = 1;     %% used for synchrosq

for rsb_target = SNR_range 
  %sigma_n = (var(s_ref) * 10^(-rsb_target/20));
  %s = s_ref + sigma_n * randn(size(s_ref));
  
  s = sigmerge(s_ref,randn(size(s_ref)),rsb_target);

  s = [z;s];                         %% 0 padding used for reconstruction
  rsb = SNR(s_ref, s(n_range0));
  fprintf('rsb_target= %f dB, rsb = %f dB\n', rsb_target, rsb);
  
  if abs(rsb - rsb_target) > 1
    error('Invalid SNR')
  end
    
  rsb = rsb_target;
   
  %fprintf(fp, '\nRecursive spectrum reconstruction : SNR=%0.02f dB\n', SNR(s(1:length(s_rec)), s_rec'));
    
 for k = k_range
  for L = L_range
  
    fprintf(1, '\n+++ Parameters: k=%d L=%0.02f snr=%0.02f M=%d +++\n', k, L, rsb, M);
    fprintf(fp, '\n+++ Parameters: k=%d L=%0.02f snr=%0.02f M=%d +++\n', k, L, rsb, M);
  
    %% 1a) Recursive STFT - index
    [tfr, nfreqs] = recursive_stft(s, k, L, mi, mf, M);
	
    %% 1b) Recursive STFT Reconstruction
    n0 = round((k-1) * L);
    [ s_hat] = stft_rec(tfr, k, L, mi,mf, M, n0);
    n_range_rec = n_range0-n0;
    s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
    fprintf(1,  '\nRecursive STFT reconstruction : RQF=%0.02f dB\n', SNR(s(n_range_rec), s_rec));
    fprintf(fp, '\nRecursive STFT reconstruction : RQF=%0.02f dB\n', SNR(s(n_range_rec), s_rec));


    %% 3a) Recursive Synchrosqueezed STFT - index3 / chemin2
    n0_opt = round((k-1) * L);
    [tfr, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0_opt);
    
    %% 3b) Reconstruction
    [ s_hat] = sstft_rec(stfr, k, L, M, n0_opt);
    n_range_rec = n_range0-n0_opt;
    s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec')))
    tmp_res = sprintf('Classical synchrosqueezing reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);
    fprintf(fp, tmp_res);

    %% n0 SNR curve
    n0_range = [n0_range0 n0_opt n0_opt+2];
    n0_curve = zeros(1, length(n0_range));
    for idx = 1:length(n0_range)
      n0r = n0_range(idx);
      [tfr, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0r);
      n_range_rec = n_range0-n0r;
      [ s_hat] = sstft_rec(stfr, k, L, M, n0r);
      s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
      n0_curve(idx) = SNR(s(n_range_rec), s_rec);
      tmp_res = sprintf('Classical synchrosqueezing reconstruction n0=%d, M=%d - RQF=%0.02f dB\n', n0r, M, n0_curve(idx));
      fprintf(1,  tmp_res);
      fprintf(fp, tmp_res);
    end
    
    %% M_range
    n0 = round((k-1) * L);
    M_curve = zeros(1, length(M_range));
    for idx = 1:length(M_range)
      Mrr = round(M_range(idx));
      [tfr, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, 1, round(Mrr/2), Mrr, n0);
      n_range_rec = n_range0-n0;
      [ s_hat] = sstft_rec(stfr, k, L, Mrr, n0);
      s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
      M_curve(idx) = SNR(s(n_range_rec), s_rec);
      tmp_res = sprintf('Classical synchrosqueezing reconstruction n0=%d, M=%d - RQF=%0.02f dB\n', n0, Mrr, M_curve(idx));
      fprintf(1,  tmp_res);
      fprintf(fp, tmp_res);
    end
    
    
    %% 3c) Recursive LM Synchrosqueezed STFT
    n0 = round((k-1) * L);
    mu_curve = zeros(1, length(mu_range));
    for idx = 1:length(mu_range)
     mu = mu_range(idx);
     [tfr, lmstfr, nfreqs, below, over] = recursive_lmsstft(s, k, L, mi, mf, M, mu, n0);
 
     %% 3d) Reconstruction
     %fprintf(1,'Signal reconstruction after synchrosqueezing\n');
     [s_hat] = sstft_rec(lmstfr, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
     mu_curve(idx) = SNR(s(n_range_rec), s_rec);
     tmp_res = sprintf('LM-Synchrosqueezing reconstruction : mu=%0.02f - RQF=%0.02f dB\n', mu,  mu_curve(idx));
     fprintf(1,  tmp_res);
     fprintf(fp, tmp_res);
     index3 = index3 + 1;
    end
   

    if plot_snr_curve
       %% n0 curve
       figure(11)
       plot(n0_range, n0_curve);
       xlabel('n0')
       ylabel('SNR (dB)');
       title(sprintf('Sy_k reconstruction quality (k=%d, L=%0.02f, M=%d, SNR=%0.02f)', k, L, M, rsb))
       saveas(gcf, sprintf('%s/n0_curve_k=%d_L=%0.02f_M=%d_SNR=%0.02f.eps', chemin3, k, L, M, rsb), 'epsc');
       
       tmp = gen_tab('$n_0$', n0_range, n0_curve);
       fprintf(1, tmp); fprintf(fp, tmp);
       
       %% M curve
       figure(12)
       plot(M_range, M_curve);
       xlabel('M')
       ylabel('SNR (dB)');
       title(sprintf('Sy_k reconstruction quality (k=%d, L=%0.02f, n0=%d, SNR=%0.02f)', k, L, n0, rsb))
       saveas(gcf, sprintf('%s/M_curve_k=%d_L=%0.02f_n0=%d_SNR=%0.02f.eps', chemin3, k, L, n0, rsb), 'epsc');
       
       tmp = gen_tab('$M$', M_range, M_curve);
       fprintf(1, tmp); fprintf(fp, tmp);
       
       %% mu curve
       figure(13)
       plot(mu_range, mu_curve);
       xlabel('\mu')
       ylabel('SNR (dB)');
       title(sprintf('Sy_k reconstruction quality (k=%d, L=%0.02f, n0=%d, M=%d, SNR=%0.02f)', k, L, n0, M, rsb))
       saveas(gcf, sprintf('%s/mu_curve_k=%d_L=%0.02f_n0=%d_M=%d_SNR=%0.02f.eps', chemin3, k, L, n0, M, rsb), 'epsc');
       tmp = gen_tab('$\\mu$', mu_range, mu_curve);
       fprintf(1, tmp); fprintf(fp, tmp);
    end
  end %% L
 end %% k
end % RSB

fclose(fp);
eps2pdf(chemin3)