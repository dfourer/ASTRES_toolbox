clear all
close all


Fs = 1000;
Ts = 1/Fs;
signal = 1;
load_signal;
s0 = s;
s0 = s0 - mean(s0);  %% zero centering

SNR_range = 45;
w0T_range = 2*pi * [2 5 7 10 15];
M_range   = [300 400 500 1000 2000];
mu_range  = [0.05 0.2 0.5 1 100];   %[0.01 0.05 0.2 0.5 1 2 5];

gamma_K     = 10^(-4);
M_default   = 500;
w0T_default = 2*pi * 7;

chemin = sprintf('figs_TSP/%d/', signal);
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

fp = fopen(sprintf('%s/results_tabMW.txt', chemin), 'w+t');

methods ={'Morlet wavelet reconstruction (Truncated Morlet formula)', 'Synchrosqueezed Morlet Wavelet reconstruction'};  %'Morlet wavelet reconstruction (Method 1)', 'S-transform reconstruction (Method 2)', 
nb_methods = length(methods);

w0      = 2*pi*100;
w_begin = 0.63;    %%% in [eps, 1]
as_range = [w0/(2*pi*w_begin) w0/(pi/Ts)];    %% [40Hz -- Fs/2Hz]


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
   T = w0T / w0;
   L = T / Ts;
   M = M_default;
   
   %[ Wx, as ] = my_MW( s, M, Ts, T, w0, as_range, gamma_K, 0)
   [ Wx, as, stfr, af, lost ] = my_sMW( s, M, Ts, T, w0, as_range, gamma_K );
   
   %% reconstruction using truncated Morlet formula
   s_hat(1,:) = 2*real(my_recMW(Wx, w0, T, as));
   
   %% synchrosqueezed Morlet WT reconstruction 
   s_hat(2,:) = real(my_recsMW(stfr, w0, T));
  
   fprintf(1, '\n');
   for i = 1:nb_methods
    w0T_tab_results(i, idx) = RQF(s.', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%.2f, f0T=%.2f, T=%.2f, L=%.2f, M=%d\tRQF=%.2f \n', methods{i}, rsb,  w0T/(2*pi), T, L, M, w0T_tab_results(i, idx));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
  %% M range
 for idx2 = 1:length(M_range)
   w0T = w0T_default;
   M   = M_range(idx2);
   
   T = w0T / w0;
   L = T / Ts;
   
   [ Wx, as, stfr, af, lost ] = my_sMW( s, M, Ts, T, w0, as_range, gamma_K );
   
   %% reconstruction using truncated Morlet formula
   s_hat(1,:) = 2*real(my_recMW(Wx, w0, T, as));
   
   %% synchrosqueezed Morlet WT reconstruction 
   s_hat(2,:) = real(my_recsMW(stfr, w0, T));
   
   fprintf(1, '\n');
   for i = 1:nb_methods
    M_tab_results(i, idx2) = RQF(s.', s_hat(i, :));
    tmp = sprintf('+ %s\t, SNR=%.2f, f0T=%.2f,T=%.2f, L=%.2f M=%d\tRQF=%.2f \n', methods{i},rsb, w0T/(2*pi), T, L, M, M_tab_results(i, idx2));
    fprintf(1,  tmp);
    fprintf(fp, tmp);
   end
 end
 
  
 for i = 1:length(methods)
  %tmp = sprintf('----  %s ---- \n', methods{i});
  %fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
  %tmp = sprintf(gen_tab('$T(s)$', w0T_range/w0, w0T_tab_results(i,:), 'a'));
  %fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  %tmp = sprintf(gen_tab('$L$', (w0T_range/w0)/Ts, w0T_tab_results(i,:), 'A'));
  %fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  tmp = sprintf(gen_tab('$f_0T$', w0T_range/(2*pi), w0T_tab_results(i,:), 'A'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);  
  tmp = sprintf(gen_tab('$M$',  M_range, M_tab_results(i,:), 'B'));
  fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
  
%   if i == 3 %% synchrosqueezing
%    tmp = sprintf(gen_tab('$\\mu$',  mu_range, mu_tab_results, 'c'));
%    fprintf(1,'%s', tmp);fprintf(fp,'%s', tmp);
%   end
  fprintf(1, '\n');
  fprintf(fp, '\n');
 end
end

fclose(fp);