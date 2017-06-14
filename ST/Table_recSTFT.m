clear all
close all

Fs = 1000;
Ts = 1/Fs;
signal = 1;
complex_signal = 0;   %% set to 1 to test reconstruction on a complex signal

load_signal;
s0 = s;
s0 = s0 - mean(s0);  %% zero centering


SNR_range = 45;
w0T_range = 2*pi * [2 5 7 10 15];
M_range   = [300 400 500 1000 2000];
mu_range  = [0.05 0.2 0.5 1 100];   %[0.01 0.05 0.2 0.5 1 2 5];

gamma_K     = 10^(-4);
M_default   = 500;
w0T_default = 7 * 2* pi;

w0 = 2*pi*100;

chemin = sprintf('figs_TSP/%d/', signal);
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

fp = fopen(sprintf('%s/results_tabGab.txt', chemin), 'w+t');

methods ={'simplified STFT reconstruction formula', 'Synch. STFT reconstruction         ', 'vertical Synch. STFT reconstruction'};
nb_methods = length(methods);


M_tab_results       = zeros(nb_methods, length(M_range));
w0T_tab_results     = zeros(nb_methods, length(w0T_range));
mu_tab_results      = zeros(1         , length(mu_range));
s_hat = zeros(nb_methods+1, length(s));

for idx0 = 1:length(SNR_range)
    
 if complex_signal == 1
   noise = hilbert(randn(size(s)));
 else
   noise = randn(size(s));
 end
 
 s = sigmerge(s0, noise, SNR_range(idx0));
 %s = s-mean(s);
 rsb = SNR(s0, s);
 
 %% w0T range
 for idx = 1:length(w0T_range)
   w0T = w0T_range(idx);
   f0T   = w0T/(2*pi);T   = w0T/w0;L   = T/Ts;
   M = M_default;
   
   %M2 = round(M/2)+1;
   %% compute the classical S-transform + synchrosqueezed versions (first- and second-order)
   [tfr, stfr]  = my_tfrsgab(s, M, L, gamma_K);
   [~, vstfr]   = my_tfrvsgab(s, M, L, gamma_K);
  
%    if complex_signal == 0    %%% if signal is real
%     %% reconstruction
%     [s_hat(1, :)] = 2*real(my_rectfrgab(tfr(1:M2,:),  L, M, m_axis(M,1:M2), M));
%     %% reconstruction SST
%     [s_hat(2, :)] = 2*real(my_rectfrsgab(stfr(1:M2,:),  L, M, m_axis(M,1:M2)));
%     %% reconstruction VSST
%     [s_hat(3, :)] = 2*real(my_rectfrsgab(vstfr(1:M2,:),  L, M, m_axis(M,1:M2)));
% 
%    else                      %% if signal is complex
%     
    %% reconstruction method 2
    [s_hat(1, :)] = my_rectfrgab(tfr,  L, M);
    %% reconstruction SST
    [s_hat(2, :)] = my_rectfrsgab(stfr,  L, M);
    %% reconstruction VSST
    [s_hat(3, :)] = my_rectfrsgab(vstfr,  L, M);
   %end
   
   fprintf(1, '\n');
   for i = 1:nb_methods
    w0T_tab_results(i, idx) = RQF(s', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%d, T=%.2f, L=%.2f, M=%d\tRQF=%.2f \n', methods{i}, rsb,  T, L, M, w0T_tab_results(i, idx));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
  %% M range
 for idx2 = 1:length(M_range)
   w0T   = w0T_default; f0T   = w0T/(2*pi);T   = w0T/w0;L   = T/Ts;
   M   = M_range(idx2);
   %M2 = round(M/2)+1;
   %% compute the classical S-transform + synchrosqueezed versions (first- and second-order)
   [tfr, stfr]  = my_tfrsgab(s, M, L, gamma_K);
   [~, vstfr]   = my_tfrvsgab(s, M, L, gamma_K);
   
%    if complex_signal == 0    %%% if signal is real
%     %% reconstruction
%     [s_hat(1, :)] = 2*real(my_rectfrgab(tfr(1:M2,:), w0T, m_axis(M,1:M2), M));
%     %% reconstruction SST
%     [s_hat(2, :)] = 2*real(my_rectfrsgab(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
%     %% reconstruction VSST
%     [s_hat(3, :)] = 2*real(my_rectfrsgab(vstfr(1:M2,:), w0T, m_axis(M,1:M2)));
% 
%    else                      %% if signal is complex
%     
   %% reconstruction method 2
   [s_hat(1, :)] = my_rectfrgab(tfr,  L, M);
   %% reconstruction SST
   [s_hat(2, :)] = my_rectfrsgab(stfr,  L, M);
   %% reconstruction VSST
   [s_hat(3, :)] = my_rectfrsgab(vstfr,  L, M);
   %end
  
   fprintf(1, '\n');
   for i = 1:nb_methods
    M_tab_results(i, idx2) = RQF(s.', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%d, T=%.2f, L=%.2f, M=%d\tRQF=%.2f \n', methods{i},rsb, T, L, M, M_tab_results(i, idx2));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
 
 for idx3 = 1:length(mu_range)
   M     = M_default;
   w0T   = w0T_default;f0T   = w0T/(2*pi);T   = w0T/w0;L   = T/Ts;
   mu    = mu_range(idx3);
   
   %% compute the LM S-transform
   [~, stfr] = my_tfrlmsgab(s, M, L, mu, gamma_K);
  
%    if complex_signal == 0    %%% if signal is real
%      [s_hat(4, :)] = 2*real(my_rectfrsgab(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
%    else
     %% reconstruction SST
     [s_hat(4, :)] = my_rectfrsgab(stfr, L, M);
%   end
   
   mu_tab_results(idx3) = RQF(s.', s_hat(4, :));
   tmp =  sprintf('+ LM-Synch STFT reconstruction\t, SNR=%d, T=%.2f, L=%.2f, M=%d, mu=%.2f\tRQF=%.2f\n ', rsb, T, L, M, mu, mu_tab_results(idx3));
   fprintf(1, tmp);
   fprintf(fp, tmp);

 end
 
 for i = 1:length(methods)
  tmp = sprintf('----  %s ---- \n', methods{i});
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  tmp = sprintf(gen_tab('$T(s)$', w0T_range/(w0), w0T_tab_results(i,:), 'a'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  tmp = sprintf(gen_tab('$L$', (w0T_range/(w0))/Ts, w0T_tab_results(i,:), 'A'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  tmp = sprintf(gen_tab('$M$',  M_range, M_tab_results(i,:), 'B'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  if i == 2 %% synchrosqueezing
   tmp = sprintf(gen_tab('$\\mu$',  mu_range, mu_tab_results, 'C'));
   fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  end
  fprintf(1, '\n');
  fprintf(fp, '\n');
 end
end

fclose(fp);