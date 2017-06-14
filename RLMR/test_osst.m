%[~, vstfr2, nfreqs, below, over, q_hat_map2] = recursive_osstft(s, k, L, mi, mf, M, n0, 2);   

[tfr, vstfr2, lost, q_hatmap, if_hatmap] = my_tfrosgab(s, M, L, 2, 2);
%[tfr, vstfr2, lost, q_hatmap, if_hatmap] = my_tfrvsgab(s, M, L, 2, 2);

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
     title(sprintf('Squared modulus of the Recurs. Obli. synchrosqueezed (CRE2t) STFT k=%d, L=%d, SNR=%d dB', k, L, rsb), 'FontSize', 14);
     %saveas(gcf, sprintf('%s/rvsstft2_%d.eps', chemin1, index), 'epsc');

          %% 3b) Reconstruction
     %[ s_hat] = sstft_rec(vstfr2, k, L, M, n0);
     [ s_hat ] = my_rectfrsgab(vstfr2, L, M).';
     
     n_range_rec = n_range0-n0;
     s_rec = real(s_hat(n_range_rec)); %% assume the signal is real
%     figure(41)
%     plot(s_rec,'k-');
%     hold on
%     plot(s, 'g-.')
%     legend('reconstruction', 'reference')
%     title(sprintf('Signal reconstruction SNR=%0.02f', SNR(s(1:length(s_rec)), s_rec)))
    tmp_res = sprintf('recursive oblique synchrosqueezing reconstruction (CRE2t) n0=(k-1)L=%0.02f, M=%d - RQF=%0.02f dB\n', n0, M, SNR(s(n_range_rec), s_rec));
    fprintf(1,  tmp_res);
    
    
    figure(223)
    plot(s(n_range_rec), 'k')
    hold on
    plot(s_rec, 'r-.')
    hold off