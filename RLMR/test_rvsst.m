% Example script for recursive Levenberg-Marquardt reassignment and recursive
% Synchrosqueezing
% 
% select the test signal setting the variable signal in 1-6
%
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

clear all
close all

sst_type = 2;  %% 1: vertical , 2:oblique


signal=1;
%% used for zero-padding for Reconstruction Quality computation
N_0pad = 100;
z      = zeros(N_0pad, 1);


%% 1: Load test signal
M = 600;
mi = 1;mf = round(M/2);

k_range   = [6 7];
L_range   = [6 8 9];
SNR_range = [40 25 8];
%n0_range0 = [8 18 27];
%mu_range = [0.3 0.8 1.3 1.8 2.3 2.8];
%M_range = [100 200 600 1000 2400];

alpha  = 0.35;           %% used for display

load_signal;
% figure(6); 
% subplot(211); plot(real(s));
% subplot(212); plot(imag(s));

n_range0 = N_0pad+(1:N);

s_ref = s;
Mr = mf-mi+1;

%q_method = 1
chemin  = sprintf('figs/', signal);

%q_method = 2;
%chemin  = sprintf('figs/KC/', signal);



chemin1 = sprintf('%s/vsstft', chemin);    %% reass
if ~exist(chemin1, 'dir')
 mkdir(chemin1);
end

t = (1:length(s))-1;
%%
index = 0;              %% used for figures indices
for rsb_target = SNR_range 

  s = sigmerge(s_ref,randn(size(s_ref)),rsb_target);
  s = [z;s];                         %% 0 padding used for reconstruction
  rsb = SNR(s_ref, s(n_range0));
  fprintf('rsb_target= %f dB, rsb = %f dB\n', rsb_target, rsb);
  
  if abs(rsb - rsb_target) > 1
    error('Invalid SNR')
  end
  rsb = rsb_target;
  
  for k = k_range
   for L = L_range
     index = index + 1;
     n0 = round((k-1) * L);
     [tfr, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0);
       

     figure(11) 
     imagesc(t, nfreqs, (abs(tfr(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     %ylabel('frequency  [Hz]', 'FontSize', 16)
     ylabel('normalized frequency', 'FontSize', 16)
     title(sprintf('Recurs. spectrogram k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/spectrogram_%d.eps', chemin1, index), 'epsc');
     
     
     figure(21) 
     imagesc(t, nfreqs, (abs(stfr(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     %ylabel('frequency  [Hz]', 'FontSize', 16)
     ylabel('normalized frequency', 'FontSize', 16)
     title(sprintf('Squared modulus of the Recurs. synchrosqueezed STFT  k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/rsstft_%d.eps', chemin1, index), 'epsc');
       
     
     %% 3b) Reconstruction
     [ s_hat] = sstft_rec(stfr, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('Classical synchrosqueezing reconstruction n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);
     
     
      %% 2a) Recursive vertical synchrosqueezing method 1
     if sst_type == 1
      [~, vstfr, nfreqs, below, over, q_hat_map] = recursive_vsstft(s, k, L, mi, mf, M, n0, 1);
     else
      [~, vstfr, nfreqs, below, over, q_hat_map] = recursive_osstft(s, k, L, mi, mf, M, n0, 1);   
     end
     figure(22)
     imagesc(t, nfreqs, (abs(vstfr(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     ylabel('normalized frequency', 'FontSize', 16)
%     %ylabel('frequency  [Hz]', 'FontSize', 16)
     title(sprintf('Squared modulus of the Recurs. Vert. synchrosqueezed (CRE1) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/rvsstft1_%d.eps', chemin1, index), 'epsc');

          %% 3b) Reconstruction
     [ s_hat] = sstft_rec(vstfr, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('Vertical synchrosqueezing reconstruction (CRE1) n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);


      %% 2a) Recursive vertical synchrosqueezing method 2
     if sst_type == 1
      [~, vstfr2, nfreqs, below, over, q_hat_map2] = recursive_vsstft(s, k, L, mi, mf, M, n0, 2);
     else
      [~, vstfr2, nfreqs, below, over, q_hat_map2] = recursive_osstft(s, k, L, mi, mf, M, n0, 2);   
     end
     figure(23)
     imagesc(t, nfreqs, (abs(vstfr2(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     ylabel('normalized frequency', 'FontSize', 16)
%     %ylabel('frequency  [Hz]', 'FontSize', 16)
     title(sprintf('Squared modulus of the Recurs. Vert. synchrosqueezed (CRE2t) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/rvsstft2_%d.eps', chemin1, index), 'epsc');

          %% 3b) Reconstruction
     [ s_hat] = sstft_rec(vstfr2, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('Vertical synchrosqueezing reconstruction (CRE2t) n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);


     %% 2b) Recursive vertical synchrosqueezing method 3 (CRE2w)
     if sst_type == 1
      [~, vstfr2, nfreqs, below, over, q_hat_map2] = recursive_vsstft(s, k, L, mi, mf, M, n0, 3);
     else
      [~, vstfr2, nfreqs, below, over, q_hat_map2] = recursive_osstft(s, k, L, mi, mf, M, n0, 3);
     end
     figure(24)
     imagesc(t, nfreqs, (abs(vstfr2(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     ylabel('normalized frequency', 'FontSize', 16)
%     %ylabel('frequency  [Hz]', 'FontSize', 16)
     title(sprintf('Squared modulus of the Recurs. Vert. synchrosqueezed (CRE2w) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/rvsstft2_%d.eps', chemin1, index), 'epsc');

          %% 3b) Reconstruction
     [ s_hat] = sstft_rec(vstfr2, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('Vertical synchrosqueezing reconstruction (CRE2w) n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);
    

    
     %% 2c) Recursive vertical synchrosqueezing method 4 (CRE2r)
     if sst_type == 1
      [~, vstfr3, nfreqs, below, over, q_hat_map3] = recursive_vsstft(s, k, L, mi, mf, M, n0, 4);
     else
      [~, vstfr3, nfreqs, below, over, q_hat_map3] = recursive_osstft(s, k, L, mi, mf, M, n0, 4);   
     end
     figure(24)
     imagesc(t, nfreqs, (abs(vstfr3(:,n_range0)).^2).^alpha);
     set(gca,'YDir','normal')
     colormap gray;
     cmap = colormap;
     cmap = flipud(cmap);
     colormap(cmap); 
     xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
     ylabel('normalized frequency', 'FontSize', 16)
%     %ylabel('frequency  [Hz]', 'FontSize', 16)
     title(sprintf('Squared modulus of the Recurs. Vert. synchrosqueezed (CRE2r) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     saveas(gcf, sprintf('%s/rvsstft2_%d.eps', chemin1, index), 'epsc');

          %% 3b) Reconstruction
     [ s_hat] = sstft_rec(vstfr3, k, L, M, n0);
     n_range_rec = n_range0-n0;
     s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('Vertical synchrosqueezing reconstruction (CRE2r) n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);
    
    
%      
%       %% qmap
%      figure(33)
%      imagesc(t, nfreqs, 10 * log10(abs(q_hat_map)));
%      %imagesc(t, nfreqs, 10 * log10(abs(q_hat_map)));
%      set(gca,'YDir','normal')
%      colormap gray;     cmap = colormap;      cmap = flipud(cmap);
%      colormap(cmap); 
%      xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
%      ylabel('normalized frequency', 'FontSize', 16)
% %     %ylabel('frequency  [Hz]', 'FontSize', 16)
%      title(sprintf('Area where q(t,\\omega) is computed (white area) is computed (CRE1) k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
%      saveas(gcf, sprintf('%s/q_hat_pbm_%d.eps', chemin1, index), 'epsc');
%      %title(sprintf('Estimated FM  k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
% 
%      figure(34)
%      imagesc(t, nfreqs, 10 * log10(abs(q_hat_map2)));
%      %imagesc(t, nfreqs, 10 * log10(abs(q_hat_map)));
%      set(gca,'YDir','normal')
%      colormap gray;     cmap = colormap;      cmap = flipud(cmap);
%      colormap(cmap); 
%      xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
%      ylabel('normalized frequency', 'FontSize', 16)
% %     %ylabel('frequency  [Hz]', 'FontSize', 16)
%      title(sprintf('Area where q(t,\\omega) is computed (white area) is computed (CRE2) k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
%      saveas(gcf, sprintf('%s/q_hat_pbm2_%d.eps', chemin1, index), 'epsc');
%      %title(sprintf('Estimated FM  k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);

     pause
   end  %% L
  end  %% k
end  % RSB
%  
eps2pdf(chemin1);