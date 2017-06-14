% Recursive versions of the Levenberg-Marquardt reassigned spectrogram 
% and of the synchrosqueezed STFT short demo
%
% Recommended requirements: 
% -TFTB (http://tftb.nongnu.org/index_fr.html)
% -ghostscript (http://ghostscript.com/download/ ps2pdf)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

clear all
close all


selected_signal = 1;  %%choose the signal to process: can be 1 or 2

if selected_signal == 1
    load fig1-observed-Livingston.txt
    t1=fig1_observed_Livingston(:,1);
    sig1=fig1_observed_Livingston(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal1';
elseif selected_signal == 2

    load fig1-observed-Hanford.txt
    t1=fig1_observed_Hanford(:,1);
    sig1=fig1_observed_Hanford(:,2);
    Te1=t1(2)-t1(1);
    Fe1=1.0/Te1;
    figs_folder = 'figs_signal2';
else
    error('Unknown signal')
end

if ~exist(figs_folder, 'dir')
 mkdir(figs_folder);
end
NbPoints=length(sig1);   %length(sig1)
alpha=0.5;             %% TFR compression factor (modify display contrast)

%% TF analysis
k_range = [2 3 4 5];                   %% recursive filter order
L_range = [50 70 120 150 220 300 400]; %70; %120;                 %% analysis window length

%% structure to store results
rqf_stft   = zeros(length(k_range), length(L_range));
rqf_sstft  = zeros(length(k_range), length(L_range));

M=5000; %mi = 0; mf = 200; %300; %1000; %
mi=floor(20*M/Fe1); mf=ceil(520*M/Fe1);   %% frequency range [20Hz - 520Hz]

n0_max = (max(k_range)-1) * max(L_range);
N_0pad = n0_max + 20; %3*250;          %% number of zeros to add before signal (used for reconstruction)
AddedZeros = zeros(N_0pad, 1);
n_range  = N_0pad+(1:length(sig1));
s = real([AddedZeros;sig1;AddedZeros]);
N = length(sig1); %NbPoints+N_0pad;         %% signal length 
t = 1:length(s)-N_0pad;

index = 0;
for i_k = 1:length(k_range)
  for i_L = 1:length(L_range)
     
    k = k_range(i_k);
    L = L_range(i_L);
    
    fprintf(1, 'k=%d L=%d\n', k, L);
    index = index + 1;
    n0 = round((k-1) * L);
    %N_0pad = n0 + 20; %3*250;          %% number of zeros to add before signal (used for reconstruction)
    %AddedZeros = zeros(N_0pad, 1);

    %% input signal
    %s = real([AddedZeros;sig1;AddedZeros]);
    %N = length(sig1); %NbPoints+N_0pad;         %% signal length     
    %t = 1:length(s)-N_0pad;

    %figure(1);
    %plot(s);

    [tfr, nfreqs] = recursive_stft(s, k, L, mi, mf, M);
    [~, rtfr, nfreqs, before, after, below, over] = recursive_rsp(s, k, L, mi, mf, M);

    %% Recursive spectrogram
    figure(2);
    imagesc(t, nfreqs, (abs(tfr(:,n_range)).^2).^alpha);
    set(gca,'YDir','normal')
    colormap gray;
    colormap(flipud(colormap)); 
    xlabel('time samples', 'FontSize', 16)
    ylabel('normalized frequency', 'FontSize', 16)
    title(sprintf('Recursive spectrogram k=%d, L=%2.2f ', k, L), 'FontSize', 14);
    saveas(gcf, sprintf('%s/rsp_%d.eps', figs_folder, index), 'epsc');

    %% Recursive reassigned spectrogram
    figure(3);
    imagesc(t, nfreqs, rtfr(:,n_range).^alpha);
    set(gca,'YDir','normal')
    colormap gray;
    colormap(flipud(colormap)); 
    xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
    ylabel('normalized frequency', 'FontSize', 16)
    title(sprintf('Recursive reassigned spectrogram k=%d, L=%2.2f ', k, L), 'FontSize', 14);
    saveas(gcf, sprintf('%s/rrsp_%d.eps', figs_folder, index), 'epsc');

    %% Recursive synchrosqueezed STFT
    %n0 = round((k-1) * L);
    [~, stfr, nfreqs, below, over] = recursive_sstft(s, k, L, mi, mf, M, n0);
    figure(4);
    imagesc(t, nfreqs, (abs(stfr(:,n_range)).^2).^alpha);
    set(gca,'YDir','normal')
    colormap gray;
    colormap(flipud(colormap)); 
    xlabel('time samples', 'FontSize', 16)  %, 'FontName', 'Times-Roman', 'FontSize', 20
    ylabel('normalized frequency', 'FontSize', 16)
    title(sprintf('Recursive synchrosqueezed STFT k=%d, L=%2.2f ', k, L), 'FontSize', 14);
    saveas(gcf, sprintf('%s/rsst_%d.eps', figs_folder, index), 'epsc');


    
    
    %% Signal reconstruction from recursive STFT
    n_range_rec = n_range;
    %n_range_rec = n_range-n0;
    s_hat = stft_rec(tfr, k, L, mi,mf, M, n0);
    s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
    rqf_stft(i_k,i_L) = RQF(s(n_range_rec), s_rec);
    fprintf(1,'Recursive STFT reconstruction : RQF=%0.02f dB\n', rqf_stft(i_k,i_L));

    figure(5);
    %plot(t(n_range_rec), s(n_range_rec), 'g-.');
    plot(s(n_range_rec), 'g-.');  %t(n_range_rec), 
    hold on
    %plot(t(n_range_rec), s_rec, 'k');
    plot(s_rec, 'k'); %t(n_range_rec), 
    legend('original signal', 'reconstructed signal')
    title('recursive STFT reconstruction')
    title(sprintf('recursive STFT reconstruction, k=%d, L=%d, M=%d, RQF=%.2f dB', k, L,M, rqf_stft(i_k,i_L)))
    hold off
    saveas(gcf, sprintf('%s/sig_stft_%d.eps', figs_folder, index), 'epsc');
    

    
    %% Signal reconstruction from recursive synchrosqueezed STFT
    s_hat = sstft_rec(stfr, k, L, M, n0);
    s_rec = 2*real(s_hat(n_range_rec)); %% assume the signal is real
    rqf_sstft(i_k,i_L) =  RQF(s(n_range_rec), s_rec);
    fprintf(1,'Recursive synchrosqueezed STFT reconstruction : RQF=%0.02f dB\n', rqf_sstft(i_k,i_L));

    figure(6);
    %plot(t(n_range_rec), s(n_range_rec), 'g-.');
    plot( s(n_range_rec), 'g-.');
    hold on
    %plot(t(n_range_rec), s_rec, 'k');
    plot( s_rec, 'k');
    legend('original signal', 'reconstructed signal')
    title(sprintf('recursive synchrosqueezed STFT reconstruction, k=%d, L=%d, M=%d, RQF=%.2f dB', k, L,M, rqf_sstft(i_k,i_L)))
    hold off
    saveas(gcf, sprintf('%s/sig_sstft_%d.eps', figs_folder, index), 'epsc');
  end
end

%% find the k which maximize RQF
[~, ik_max] = max(mean(rqf_stft.'));

%% plot 

fprintf(1,'Best results for recursive STFT reconstruction with k=%d\n', k_range(ik_max))
tmp = gen_tab('$L$', L_range, rqf_stft(ik_max,:));
fprintf(1, tmp);

%% find the k which maximize RQF
[~, ik_max2] = max(mean(rqf_sstft.'));

%% plot 
fprintf(1,'Best results for recursive synchrosqueezed STFT reconstruction with k=%d\n', k_range(ik_max2))
tmp = gen_tab('$L$', L_range, rqf_sstft(ik_max2,:));
fprintf(1, tmp);


eps2pdf(figs_folder);


