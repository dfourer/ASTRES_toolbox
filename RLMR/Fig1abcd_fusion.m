%
%  Generate Fig 1 (a)-(d) of Ref
%  Performance curves of Chirp Rate (CR) and Instantaneous Frequency (IF)
%  estimators
%  Merge results with amplitude modulation
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 27-06-2016
% Ref: [D. Fourer, F. Auger, K.Czarnecki, S. Meignen and P. Flandrin, Chirp rate and instantaneous frequency estimation. Proc. IEEE ICASSP 2017]

clear all
close all

chemin = 'figs';
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

use_recursive_stft = true;

%                  IF   CRE
processed_cases = [1    1;...   %% IF(1) - CRE1
                   2    2;...   %% IF(2) - CRE2t
                   2    3;...   %% IF(2) - CRE2w
                   2    4;...   %% IF(2) - CRE2r
                   3    1];     %% IF(3) - x
nb_cases        = size(processed_cases,1);


%% files to use
if use_recursive_stft
 matfilename1 = 'Figabcd_recursive.mat';
else
 matfilename1 = 'Figabcd_gaussian.mat';
end

matfilename2 = sprintf('Tx_%s', matfilename1);

if ~exist(matfilename1, 'file') ||  ~exist(matfilename2, 'file')
  error('Please launch Fig1abcd_comparison first with use_recursive_stft=%d', use_recursive_stft);
end

load(matfilename1);
xxx = load(matfilename2, 'q_results', 'if_results');
q_results2 = xxx.q_results;
if_results2 = xxx.if_results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 1a MSE(CRE) as a function of frequency

%q_methods = {'CRE1' 'CRE2t' 'CRE2\omega' 'CRE2r'}; % };
q_methods =  {'$\hat{\alpha}_x^{K1} (L_x=+\infty)$' '$\hat{\alpha}_x^{t2} (L_x=+\infty)$' '$\hat{\alpha}_x^{\omega 2} (L_x=+\infty)$' '$\hat{\alpha}_x^{r2} (L_x=+\infty)$'}; % };
q_methods2 = {'$\hat{\alpha}_x^{K1} (L_x=200)$'     '$\hat{\alpha}_x^{t2} (L_x=200)$'     '$\hat{\alpha}_x^{\omega 2} (L_x=200)$'     '$\hat{\alpha}_x^{r2} (L_x=200)$'}; % };
cols = {'k-s', 'r-^', 'b-v', 'g-x','m-o', 'r-v', 'b-o'};
cols2 = {'k-.s', 'r-.^', 'b-.v', 'g-.x','m-.o', 'r-.v', 'b-.o'};

for idx_snr = 1:length(snr_range)
 figure(1)
 index = 1;
 for c = 1:4
  if_method = processed_cases(c,1);
  q_method  = processed_cases(c,2);
  
  mse = mean( reshape( q_results(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  %mse = mean( reshape( q_results(idx_snr, q_method, 1,1, :), nb_iter, nb_freq),1);
  
  %plot(n_freq(nfr), mse(nfr), cols{q_method});
  %xlabel('normalized frequency', 'FontSize', 16);
  %ylabel('MSE', 'FontSize', 16)
  
  l(index)   = plot(n_freq, 10*log10(mse+eps), cols{q_method});
  leg{index} = q_methods{c};
  xlabel('normalized frequency \lambda', 'FontSize', 16);
  ylabel('MSE [dB]', 'FontSize', 16)
  
  hold on
  index = index + 1;
 end
 
 for c = 1:4
  if_method = processed_cases(c,1);
  q_method  = processed_cases(c,2);
  
  mse = mean( reshape( q_results2(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  %mse = mean( reshape( q_results(idx_snr, q_method, 1,1, :), nb_iter, nb_freq),1);
  
  %plot(n_freq(nfr), mse(nfr), cols{q_method});
  %xlabel('normalized frequency', 'FontSize', 16);
  %ylabel('MSE', 'FontSize', 16)
  
  l(index)   = plot(n_freq, 10*log10(mse+eps), cols2{q_method});
  leg{index} = q_methods2{c};
  xlabel('normalized frequency \lambda', 'FontSize', 16);
  ylabel('MSE [dB]', 'FontSize', 16)
  
  hold on
  index = index + 1;
 end
 
 ylim([-120 0]); grid;
 title(sprintf('CRE, SNR=%.2f dB', snr_range(idx_snr)), 'FontSize', 16)
 legend(l, leg, 'Location', 'NorthWest', 'Interpreter','Latex'); %'Best'
 
 saveas(gcf, sprintf('%s/fig1a_mse_cre_%d.eps', chemin,idx_snr), 'epsc');
 hold off
end
eps2pdf(chemin);



freq_range = {[1:length(n_freq)], [11], [26]};

figure(2)
for idx = 1:length(freq_range)

  index = 1;
  mse_snr = zeros(1, length(snr_range));
  for c = 1:4
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    for id_snr = 1:length(snr_range)
      tmp = reshape(q_results(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, nb_iter*length(freq_range{idx}));   %% squared error
            
      mse_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp);   %sum(sum(q_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq);
      %mse_snr(id_snr) = mean(mean(reshape(q_results(id_snr, :, q_method, 1, :).^2, nb_iter, nb_freq),1));
    end
    l(index)    = plot(snr_range, 10 * log10(eps+mse_snr), cols{c});
    leg{index}  = q_methods{c};
    hold on
    index = index +1;
  end
  for c = 1:4
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    for id_snr = 1:length(snr_range)
      tmp = reshape(q_results2(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, nb_iter*length(freq_range{idx}));   %% squared error
            
      mse_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp);   %sum(sum(q_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq);
      %mse_snr(id_snr) = mean(mean(reshape(q_results(id_snr, :, q_method, 1, :).^2, nb_iter, nb_freq),1));
    end
    l(index)    = plot(snr_range, 10 * log10(eps+mse_snr), cols2{c});
    leg{index}  = q_methods2{c};
    index = index +1;
  end
  
  if length(freq_range{idx}) == 1
    title(sprintf('CRE, $\\lambda=%.2f$',n_freq(freq_range{idx})), 'Interpreter','Latex','FontSize', 16);
  else
    title(sprintf('CRE, $\\lambda\\in[%.2f,%.2f]$',n_freq(freq_range{idx}(1)),n_freq(freq_range{idx}(end))), 'Interpreter','Latex','FontSize', 16); 
  end
  ylim([-120 0]); grid;
  legend(l, leg, 'Location', 'SouthWest', 'Interpreter','Latex');
  xlabel('SNR [dB]','FontSize', 16);
  ylabel('MSE [dB]','FontSize', 16);
  saveas(gcf, sprintf('%s/fig1b_mse_cre_snr_%d.eps', chemin, idx), 'epsc');
  hold off
end
eps2pdf(chemin);



%% Fig 1c MSE(IF) as a function of frequency
%if_methods = {'(1)' '(2)' 'classic'};
%q_methods = {'CRE1' 'CRE2t' 'CRE2\omega' 'CRE2r'}; % };

ifq_methods  = {'$\hat{\dot{\phi}}_x^{K1} (L_x=+\infty)$' '$\hat{\dot{\phi}}_x^{t2} (L_x=+\infty)$' '$\hat{\dot{\phi}}_x^{\omega 2} (L_x=+\infty)$' '$\hat{\dot{\phi}}_x^{r2} (L_x=+\infty)$' '$\hat{\omega} (L_x=+\infty)$'};
ifq_methods2 = {'$\hat{\dot{\phi}}_x^{K1} (L_x=200)$'     '$\hat{\dot{\phi}}_x^{t2} (L_x=200)$'     '$\hat{\dot{\phi}}_x^{\omega 2} (L_x=200)$'     '$\hat{\dot{\phi}}_x^{r2} (L_x=200)$'     '$\hat{\omega} (L_x=200)$'};
%cols = {'k-', 'k-^', 'r-.', 'r-x','g-o', 'g-v', 'c-v'};

for idx_snr = 1:length(snr_range)
  figure(3)
  index = 1;
  for c = 1:nb_cases
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    %fprintf(1, '%s \n', if_methods{if_method});

    mse = mean( reshape( if_results(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  
    l(index)   = plot(n_freq, 10*log10(mse+eps), cols{c});
    leg{index} = sprintf('%s', ifq_methods{c});
    hold on
    index = index + 1;
  end
  for c = 1:nb_cases
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    %fprintf(1, '%s \n', if_methods{if_method});

    mse = mean( reshape( if_results2(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  
    l(index)   = plot(n_freq, 10*log10(mse+eps), cols2{c});
    leg{index} = sprintf('%s', ifq_methods2{c});
    hold on
    index = index + 1;
  end;
  
 grid; 
 ylim([-120 0]);
 title(sprintf('IFE, SNR=%.2f dB', snr_range(idx_snr)), 'FontSize', 16)
 legend(l, leg, 'Location', 'NorthWest', 'Interpreter', 'Latex');
 xlabel('normalized frequency \lambda','FontSize', 16);
 ylabel('MSE [dB]','FontSize', 16)
 saveas(gcf, sprintf('%s/fig1c_mse_if_%d.eps', chemin, idx_snr), 'epsc');
 hold off
end
eps2pdf(chemin);


%% Fig 1d MSE(IF) as a function of SNR
figure(4)

mseif_snr = zeros(1, length(snr_range));
for idx = 1:length(freq_range)
 index = 1;
 for c = 1:nb_cases
   if_method = processed_cases(c,1);
   q_method  = processed_cases(c,2);
   
   for id_snr = 1:length(snr_range)
     tmp = reshape(if_results(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, (nb_iter * length(freq_range{idx})));
     
     mseif_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp); %
     %sum(sum(if_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq); 
   end
   if length(freq_range{idx}) == 1
    title(sprintf('IFE, $\\lambda=%.2f$',n_freq(freq_range{idx})), 'Interpreter','Latex','FontSize', 16);
   else
    title(sprintf('IFE, $\\lambda\\in[%.2f,%.2f]$',n_freq(freq_range{idx}(1)),n_freq(freq_range{idx}(end))), 'Interpreter','Latex','FontSize', 16); 
   end
   l(index)   = plot(snr_range, 10*log10(mseif_snr+eps), cols{c});
   leg{index} = ifq_methods{c};
   hold on
   index = index + 1;
 end
 for c = 1:nb_cases
   if_method = processed_cases(c,1);
   q_method  = processed_cases(c,2);
   
   for id_snr = 1:length(snr_range)
     tmp = reshape(if_results2(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, (nb_iter * length(freq_range{idx})));
     
     mseif_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp); %
     %sum(sum(if_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq); 
   end
   if length(freq_range{idx}) == 1
    title(sprintf('IFE, $\\lambda=%.2f$',n_freq(freq_range{idx})), 'Interpreter','Latex','FontSize', 16);
   else
    title(sprintf('IFE, $\\lambda\\in[%.2f,%.2f]$',n_freq(freq_range{idx}(1)),n_freq(freq_range{idx}(end))), 'Interpreter','Latex','FontSize', 16); 
   end
   l(index)   = plot(snr_range, 10*log10(mseif_snr+eps), cols2{c});
   leg{index} = ifq_methods2{c};
   index = index + 1;
 end
 ylim([-120 0]); grid;
 legend(l, leg,  'Interpreter','Latex');
 xlabel('SNR [dB]','FontSize', 16);
 ylabel('MSE [dB]','FontSize', 16);
 legend(l, leg, 'Location', 'SouthWest', 'Interpreter', 'Latex');

 %legend(l, leg, 'Location', 'SouthWest');  %NorthEast
 saveas(gcf, sprintf('%s/fig1d_mse_if_snr_%d.eps', chemin, idx), 'epsc');
 hold off
 
end

eps2pdf(chemin);




