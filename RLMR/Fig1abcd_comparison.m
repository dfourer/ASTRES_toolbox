%
%  Generate Fig 1 (a)-(d) of Ref
%  Performance curves of Chirp Rate (CR) and Instantaneous Frequency (IF)
%  estimators
%
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

recompute_all = false;     % force to recompute all neither if matfilename exist
use_recursive_stft = true; % use the recursive implementation or not
amp_modulation = true;    % consider amplitude modulation or not

if use_recursive_stft
 matfilename = 'Figabcd_recursive.mat';
else
 matfilename = 'Figabcd_gaussian.mat';
end

if amp_modulation
  matfilename = sprintf('Tx_%s', matfilename);
end


%                  IF   CRE
processed_cases = [1    1;...   %% IF(1) - CRE1
                   2    2;...   %% IF(2) - CRE2t
                   2    3;...   %% IF(2) - CRE2w
                   2    4;...   %% IF(2) - CRE2r
                   3    1];     %% IF(3) - x
nb_cases        = size(processed_cases,1);

if recompute_all || ~exist(matfilename, 'file')

M = 500;
Mh = round(M/2);

L = 6;  %6; % [5 7 10];
gamma_K = 1e-4;
q_threshold = 1e-4;

nb_qmethods  = 4;
nb_ifmethods = 3;

nb_iter = 50; %100;   %% 10 %% 20
snr_range = 80:(-5):(-10);  % [%[80 70 60 50 40 30 20 10 5 0 -5 -10];

idx_snr = 2;

if use_recursive_stft
  k = 7;
  n0 = (k-1)*L;
end;

mi = 0;
mf = Mh;
  
%% generate signal
t0 = 250; N = 500;
alpha = 2*pi*0.36/N;
t = (0:N-1);
A   = 1;
if amp_modulation
 T_x = 200;   
else
 T_x = inf; % 200; %200;   %% set T_x to inf for constant amplitude
end
Ax    = A * exp(-(t-t0).^2 / (2*T_x^2));
phi_x =  alpha * t.^2/2;  %2*pi * 0 *t +

s = Ax .* exp(1j * phi_x);   % + real(Ax .* exp(1j .* 2*pi * 0.4 * t));
s = s(:);

%% center frequency  
lambda_0 = t0 * alpha/(2*pi);    %% d phi_x / dt
m_0 = floor(lambda_0 * M)+1;     %% discrete frequency bin (Matlab index)

if use_recursive_stft
 %dM = 25;
 m_rg = ((m_0-25):(m_0+10))-mi;
else
 dM = 36;
 m_rg = ((m_0-dM):(m_0+dM))-mi;
end

n_freq = (m_rg-1)/M;      %% computed frequency axis
nb_freq = length(n_freq); %% number of bins 2*dM+1;

q_results  = zeros(length(snr_range), nb_iter, nb_qmethods, nb_ifmethods, nb_freq);
if_results = zeros(length(snr_range), nb_iter, nb_qmethods, nb_ifmethods, nb_freq);

index = 0;
for rsb = snr_range
 index = index + 1;
 for id = 1:nb_iter
   tic
   fprintf(1,'Iteration %d / %d, SNR=%.2f \n', id, nb_iter, rsb);
   x = sigmerge(s, hilbert(randn(size(s))), rsb); %% add noise
   
   %% process each case
   qr_tmp  = zeros(nb_cases, length(m_rg));
   ifr_tmp = zeros(nb_cases, length(m_rg));
   parfor c = 1:nb_cases
   
       if_m = processed_cases(c,1);
       q_m  = processed_cases(c,2);
       
%        
%    for q_method = 1:nb_qmethods
%      %m_rg = m_rg(m_rg>=1);
%      %m_rg = m_rg(m_rg<=M);
%      
%      parfor if_method = 1:nb_ifmethods
%       
%       fprintf(1,'Computing q_method %d/%d, if_method %d/%d \n', q_method, nb_qmethods, if_method, nb_ifmethods);
%                
      if use_recursive_stft
       [tfr, stfr, nfreqs, below, over, q_hatmap, if_hatmap] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_m, if_m, q_threshold);
      else
       [tfr, stfr, lost, q_hatmap, if_hatmap] = my_tfrvsgab(x, M, L, q_m, if_m, gamma_K, q_threshold);
      end
%     
% %       %% debug (uncomment for unitary test)
% %       contraste = 0.35;
% %       figure(1)
% %       plot(n_freq, (abs(stfr(m_rg, t0+1)).^2).^contraste, 'k-')
% %       hold on
% %       plot(n_freq, (abs(tfr(m_rg, t0+1)).^2).^contraste, 'g-.')
% %       plot((m_0-1)/M, (abs(stfr(m_0, t0+1)).^2).^contraste, 'ro')
% %       xlabel('\lambda')
% %       ylabel('amplitude')
% %       legend('squ. mod. of recursiv. vert. synch. STFT', 'recursive spectrogram', 'true IF', 'Location', 'Best')
% %       title(sprintf('VSST CRE%d IF%d', q_method, if_method));
% %       if use_recursive_stft, prefix = 'recursive'; else, prefix=''; end
% %       saveas(gcf, sprintf('test_unitaires/%s_VSSTFT_CRE%d_IF%d.eps', prefix, q_method, if_method));
% %       hold off
% %       %% 
% %       figure(2)
% %       plot(n_freq, alpha * ones(1, length(m_rg)), 'r-.')
% %       hold on
% %       plot(n_freq, imag(q_hatmap(m_rg,t0+1)))
% %       title(sprintf('CRE%d', q_method));
% %       xlabel('\lambda')
% %       ylabel('CRE [cycle/sample^2]')
% %       saveas(gcf, sprintf('test_unitaires/%s_CRE%d.eps', prefix, q_method));
% %       hold off
% %       
% %       figure(3)
% %       plot(n_freq, lambda_0 * ones(1, length(m_rg)), 'r-.')
% %       hold on
% %       plot(n_freq, if_hatmap(m_rg,t0+1))
% %       xlabel('\lambda')
% %       ylabel('IF [cycle/sample]')
% %       title(sprintf('IF%d', if_method));
% %       saveas(gcf, sprintf('test_unitaires/%s_IF%d.eps', prefix, if_method));
% %       hold off
% %       
% %       if_hatmap(m_0,t0+1)-lambda_0
% %       pause
% %       %% end debug
% 
        qr_tmp(c,:)  = imag(q_hatmap(m_rg,t0+1))-alpha;
        ifr_tmp(c,:) = if_hatmap(m_rg,t0+1) -lambda_0;
%      end
   end
   
   for c = 1:nb_cases
    if_m = processed_cases(c,1);
    q_m  = processed_cases(c,2);  
    q_results(index , id, q_m, if_m, :) = qr_tmp(c,:);
    if_results(index, id, q_m, if_m, :) = ifr_tmp(c,:);
   end
   
   
   toc
 end
end
  save(matfilename);
else
  load(matfilename);
end


%% Fig 1a MSE(CRE) as a function of frequency

%q_methods = {'CRE1' 'CRE2t' 'CRE2\omega' 'CRE2r'}; % };
q_methods = {'$\hat{\alpha}_x^{K1}$' '$\hat{\alpha}_x^{t2}$' '$\hat{\alpha}_x^{\omega 2}$' '$\hat{\alpha}_x^{r2}$'}; % };
cols = {'k-s', 'r-^', 'b-v', 'g-x','y-o', 'r-v', 'b-o'};

for idx_snr = 1:length(snr_range)
 figure(1)
 for c = 1:4
  if_method = processed_cases(c,1);
  q_method  = processed_cases(c,2);
  
  %fprintf(1, '%d - %d\n', if_method, q_method);
  mse = mean( reshape( q_results(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  %mse = mean( reshape( q_results(idx_snr, q_method, 1,1, :), nb_iter, nb_freq),1);
  
  %plot(n_freq(nfr), mse(nfr), cols{q_method});
  %xlabel('normalized frequency', 'FontSize', 16);
  %ylabel('MSE', 'FontSize', 16)
  
  plot(n_freq, 10*log10(mse+eps), cols{q_method});
  xlabel('normalized frequency \lambda', 'FontSize', 16);
  ylabel('MSE [dB]', 'FontSize', 16)
  
  hold on
 end
 ylim([-110 0])
 title(sprintf('CRE, SNR=%.2f dB', snr_range(idx_snr)), 'FontSize', 16)
 legend(q_methods, 'Location', 'NorthWest', 'Interpreter','Latex'); %'Best'
 
 saveas(gcf, sprintf('%s/fig1a_mse_cre_%d.eps', chemin,idx_snr), 'epsc');
 hold off
end
eps2pdf(chemin);


%% Fig 1b MSE(CRE) as a function of SNR
%1: \lambda=0.15   2: \lambda:0.18
freq_range = {[1:length(n_freq)], [11], [26]};

figure(2)
for idx = 1:length(freq_range)

  mse_snr = zeros(1, length(snr_range));

  for c = 1:4
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    for id_snr = 1:length(snr_range)
      tmp = reshape(q_results(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, nb_iter*length(freq_range{idx}));   %% squared error
            
      mse_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp);   %sum(sum(q_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq);
      %mse_snr(id_snr) = mean(mean(reshape(q_results(id_snr, :, q_method, 1, :).^2, nb_iter, nb_freq),1));
    end
    plot(snr_range, 10 * log10(eps+mse_snr), cols{q_method});
    hold on
  end
  if length(freq_range{idx}) == 1
    title(sprintf('$\\lambda=%.2f$',n_freq(freq_range{idx})), 'Interpreter','Latex','FontSize', 16);
  else
    title(sprintf('$\\lambda\\in[%.2f,%.2f]$',n_freq(freq_range{idx}(1)),n_freq(freq_range{idx}(end))), 'Interpreter','Latex','FontSize', 16); 
  end
  legend(q_methods, 'Interpreter','Latex','Location', 'SouthWest');
  xlabel('SNR','FontSize', 16);
  ylabel('MSE [dB]','FontSize', 16);
  saveas(gcf, sprintf('%s/fig1b_mse_cre_snr_%d.eps', chemin, idx), 'epsc');
  hold off
end
eps2pdf(chemin);


%%%%%%%%%%%%%%%%%%%%%%%  Frequency estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fig 1c MSE(IF) as a function of frequency
%if_methods = {'(1)' '(2)' 'classic'};
%q_methods = {'CRE1' 'CRE2t' 'CRE2\omega' 'CRE2r'}; % };

ifq_methods = {'$\hat{\dot{\phi}}_x^{K1}$' '$\hat{\dot{\phi}}_x^{t2}$' '$\hat{\dot{\phi}}_x^{\omega 2}$' '$\hat{\dot{\phi}}_x^{r2}$' '$\hat{\omega}$'};

%cols = {'k-', 'k-^', 'r-.', 'r-x','g-o', 'g-v', 'c-v'};
l = [];
leg = {};

for idx_snr = 1:length(snr_range)
  figure(3)
  for c = 1:nb_cases
    if_method = processed_cases(c,1);
    q_method  = processed_cases(c,2);
    %fprintf(1, '%s \n', if_methods{if_method});

    mse = mean( reshape( if_results(idx_snr, :, q_method, if_method, :).^2, nb_iter, nb_freq),1);
  
    l(c)   = plot(n_freq, 10*log10(mse+eps), cols{c});
    leg{c} = sprintf('%s', ifq_methods{c});
    hold on
  end
 ylim([-90 10])
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
l = [];
leg = {};
mseif_snr = zeros(1, length(snr_range));
for idx = 1:length(freq_range)
 for c = 1:nb_cases
   if_method = processed_cases(c,1);
   q_method  = processed_cases(c,2);
   
   for id_snr = 1:length(snr_range)
     tmp = reshape(if_results(id_snr, :, q_method, if_method, freq_range{idx}).^2, 1, (nb_iter * length(freq_range{idx})));
     
     mseif_snr(id_snr) = robust_mean(tmp); %robust_mean(tmp); %
     %sum(sum(if_results(id_snr, :, q_method, if_method, :).^2)) / (nb_iter * nb_freq); 
   end
   if length(freq_range{idx}) == 1
    title(sprintf('$\\lambda=%.2f$',n_freq(freq_range{idx})), 'Interpreter','Latex','FontSize', 16);
   else
    title(sprintf('$\\lambda\\in[%.2f,%.2f]$',n_freq(freq_range{idx}(1)),n_freq(freq_range{idx}(end))), 'Interpreter','Latex','FontSize', 16); 
   end
   l(c)   = plot(snr_range, 10*log10(mseif_snr+eps), cols{c});
   leg{c} = ifq_methods{c};
   hold on
 end
 legend(l, leg,  'Interpreter','Latex');
 xlabel('SNR','FontSize', 16);
 ylabel('MSE [dB]','FontSize', 16);
 legend(l, leg, 'Location', 'SouthWest', 'Interpreter', 'Latex');

 %legend(l, leg, 'Location', 'SouthWest');  %NorthEast
 saveas(gcf, sprintf('%s/fig1d_mse_if_snr_%d.eps', chemin, idx), 'epsc');
 hold off
end

eps2pdf(chemin);

