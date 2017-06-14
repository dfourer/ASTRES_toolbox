clear all
close all

Fs = 1000;
signal = 1;
complex_signal = 0;   %% set to 1 to test reconstruction on a complex signal

fast_processing = 0;  %% if 1, computation time is much longer, if 0, please ignore results for ST reconstruction method 1

load_signal;
s0 = s;
s0 = s0 - mean(s0);  %% zero centering


SNR_range = 45;
w0T_range = 2*pi * [2 5 7 10 15];
M_range   = [300 400 500 1000 2000];   %100 200
mu_range  = [0.05 0.2 0.5 1 100];   %[0.01 0.05 0.2 0.5 1 2 5];

gamma_K     = 10^(-4);
M_default   = 500;
w0T_default = 7 * 2* pi;

chemin = sprintf('figs_TSP/%d/', signal);
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

fp = fopen(sprintf('%s/results_tab.txt', chemin), 'w+t');

methods ={'S-transform reconstruction (Method 1)', 'S-transform reconstruction (Method 2)', 'Synch. ST reconstruction         ', 'Vertical Synch. ST reconstruction'};
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
 s = s-mean(s);
 rsb = SNR(s0, s);
 
 
 
 %% w0T range
 for idx = 1:length(w0T_range)
   w0T = w0T_range(idx);   
   f0T = w0T/(2*pi);
   M = M_default;
   %M2 = round(M/2)+1;
   %% compute the classical S-transform + synchrosqueezed versions (first- and second-order)
   
   %z_pad  = round(w0T * M * sqrt(2*log(1/gamma_K)) / pi);
   %s0p = [zeros(z_pad,1) ; s ; zeros(z_pad, 1)];
   %i_s = z_pad+1:(length(s0p)-z_pad);
 
   if fast_processing == 1
    [tfr0p] = [zeros(M,1) tfrst(s, M, w0T, gamma_K) zeros(M,1)];
   else
    [tfr0p] = tfrstzp(s, M, w0T, gamma_K);
   end
   
   %[tfr0p]      =  tfrst(s0p, M, w0T, gamma_K);
   [tfr, stfr]  =  tfrsst(s, M, w0T, gamma_K);
   [~,   vstfr] = tfrvsst(s, M, w0T, gamma_K);
  
%    if complex_signal == 0    %%% if signal is real
%     %% reconstruction method 1
%     stmp = 2*real(rectfrst(tfr0p(1:M2,:), m_axis(M,1:M2), M));
%     %% reconstruction method 2
%     [s_hat(2, :)] = 2*real(rectfrst2(tfr(1:M2,:), w0T, m_axis(M,1:M2), M));
%     %% reconstruction SST
%     [s_hat(3, :)] = 2*real(rectfrsst(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
%     %% reconstruction VSST
%     [s_hat(4, :)] = 2*real(rectfrsst(vstfr(1:M2,:), w0T, m_axis(M,1:M2)));
% 
%    else                      %% if signal is complex
    
    %% reconstruction method 1
    stmp = rectfrst(tfr0p);
    %% reconstruction method 2
    [s_hat(2, :)] = rectfrst2(tfr, w0T);
    %% reconstruction SST
    [s_hat(3, :)] = rectfrsst(stfr, w0T);
    %% reconstruction VSST
    [s_hat(4, :)] = rectfrsst(vstfr, w0T);
%   end
   [s_hat(1, :)] = stmp(2:(end-1));
%    plot(real(s_hat(1, :)))
%    hold on
%    plot(real(s),'r-.')
%    pause
   fprintf(1, '\n');
   for i = 1:nb_methods
    w0T_tab_results(i, idx) = RQF(s.', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%d, w0T=%.2f, f0T=%.2f, M=%d\tRQF=%.2f \n', methods{i}, rsb,  w0T, f0T, M, w0T_tab_results(i, idx));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
  %% M range
 for idx2 = 1:length(M_range)
   w0T = w0T_default;
   f0T = w0T/(2*pi);
   M   = M_range(idx2);
   %M2 = round(M/2)+1;
   
   %% compute the classical S-transform + synchrosqueezed versions (first- and second-order)
   if fast_processing == 1
    [tfr0p] = [zeros(M,1) tfrst(s, M, w0T, gamma_K) zeros(M,1)];
   else
    [tfr0p] = tfrstzp(s, M, w0T, gamma_K);
   end
   [tfr, stfr]  =  tfrsst(s, M, w0T, gamma_K);
   [~,   vstfr] = tfrvsst(s, M, w0T, gamma_K);
  
%    if complex_signal == 0    %%% if signal is real
%     %% reconstruction method 1
%     stmp = 2*real(rectfrst(tfr0p(1:M2,:), m_axis(M,1:M2), M));
%     %% reconstruction method 2
%     [s_hat(2, :)] = 2*real(rectfrst2(tfr(1:M2,:), w0T, m_axis(M,1:M2), M));
%     %% reconstruction SST
%     [s_hat(3, :)] = 2*real(rectfrsst(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
%     %% reconstruction VSST
%     [s_hat(4, :)] = 2*real(rectfrsst(vstfr(1:M2,:), w0T, m_axis(M,1:M2)));
% 
%    else                      %% if signal is complex
    
    %% reconstruction method 1
    stmp = rectfrst(tfr0p);
    %% reconstruction method 2
    [s_hat(2, :)] = rectfrst2(tfr, w0T);
    %% reconstruction SST
    [s_hat(3, :)] = rectfrsst(stfr, w0T);
    %% reconstruction VSST
    [s_hat(4, :)] = rectfrsst(vstfr, w0T);
%   end
   [s_hat(1, :)] = stmp(2:(end-1));
   fprintf(1, '\n');
   
   for i = 1:nb_methods
    M_tab_results(i, idx2) = RQF(s.', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%d, w0T=%.2f, f0T=%.2f, M=%d\tRQF=%.2f \n', methods{i},rsb,  w0T,f0T, M, M_tab_results(i, idx2));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
 
 for idx3 = 1:length(mu_range)
   M     = M_default;
   %M2 = round(M/2)+1;
   
   w0T   = w0T_default;
   f0T   = w0T/(2*pi);
   mu    = mu_range(idx3);
   
   %% compute the LM S-transform
   [tfr, stfr]  = tfrlmsst(s, mu,M, w0T, gamma_K);
  
   %if complex_signal == 0    %%% if signal is real
   %  [s_hat(5, :)] = 2*real(rectfrsst(stfr(1:M2,:), w0T, m_axis(M,1:M2)));
   %else
     %% reconstruction SST
     [s_hat(5, :)] = rectfrsst(stfr, w0T);
   %end
   
   mu_tab_results(idx3) = RQF(s.', s_hat(5, :));
   tmp =  sprintf('+ LM-Synch ST reconstruction\t, SNR=%d, w0T=%.2f, f0T=%.2f, M=%d, mu=%.2f\tRQF=%.2f\n ', rsb, w0T, f0T, M, mu, mu_tab_results(idx3));
   fprintf(1, tmp);
   fprintf(fp, tmp);

 end
 
 for i = 1:length(methods)
  tmp = sprintf('----  %s ---- \n', methods{i});
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  tmp = sprintf(gen_tab('$f_{0}T$', w0T_range/(2*pi), w0T_tab_results(i,:), 'A'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  tmp = sprintf(gen_tab('$M$',  M_range, M_tab_results(i,:), 'B'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  if i == 3 %% synchrosqueezing
   tmp = sprintf(gen_tab('$\\mu$',  mu_range, mu_tab_results, 'C'));
   fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  end
  fprintf(1, '\n');
  fprintf(fp, '\n');
 end
end

fclose(fp);