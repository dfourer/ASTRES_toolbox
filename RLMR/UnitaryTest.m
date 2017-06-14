chemin = 'test_unitaires';
if ~exist(chemin, 'dir')
mkdir(chemin);
end

use_recursive_stft = true; 

M = 500; Mh = round(M/2);
L = 6;  %6; % [5 7 10];
q_threshold = 1e-4;

k = 7; %6;
n0 = (k-1)*L;

mi = 0;
mf = Mh;

%% generate signal
N = 500;
alpha = 2*pi*0.36/N;
t = (0:N-1); t0 = 250;

A     = 1;
T_x   = inf; % 200; %200;   %% set T_x to inf for constant amplitude
Ax    = A * exp(-(t-t0).^2 / (2*T_x^2));
phi_x =  alpha * t.^2/2;  %2*pi * 0 *t +

s = Ax .* exp(1j * phi_x);   % + real(Ax .* exp(1j .* 2*pi * 0.4 * t));
s = s(:);

lambda_0 = t0 * alpha/(2*pi);    %% d phi_x / dt
m_0 = floor(lambda_0 * M)+1;     %% discrete frequency bin (Matlab index)
m_rg = (m_0-25):(m_0+10)-mi;

rsb = 100;
x = sigmerge(s, hilbert(randn(size(s))), rsb); %% add noise  

q_method = 2; 
q_methods = {'CRE1' 'CRE2t' 'CRE2\omega' 'CRE2r'}; % };
if_method = 2;

[tfr, stfr, nfreqs, below, over, q_hatmap, if_hatmap] = recursive_vsstft(x, k, L, mi, mf, M, n0, q_method, if_method, q_threshold);       
n_freq  = nfreqs(m_rg);   %% computed frequency axis
nb_freq = length(n_freq); %% number of bins 2*dM+1;

contraste = 0.35;
figure(1)
plot(n_freq, (abs(stfr(m_rg, t0+1)).^2).^contraste, 'k-',...
     n_freq, (abs(tfr(m_rg, t0+1)).^2).^contraste, 'g-.',...
     (m_0-1)/M, (abs(stfr(m_0, t0+1)).^2).^contraste, 'ro');
xlabel('\lambda'); ylabel('amplitude');
legend('squ. mod. of recursiv. vert. synch. STFT', 'recursive spectrogram', 'true IF', 'Location', 'Best')
title(sprintf('VSST CRE%d IF%d', q_method, if_method));
if use_recursive_stft, prefix = 'recursive'; else, prefix=''; end
saveas(gcf, sprintf('%s/%s_VSSTFT_CRE%d_IF%d.eps', chemin, prefix, q_method, if_method));

%% 
figure(2); plot(n_freq, imag(q_hatmap(m_rg,t0+1))-alpha); grid;
title(sprintf('CRE%d', q_method));
xlabel('\lambda'); ylabel('CRE estimation error [cycle/sample^2]')
saveas(gcf, sprintf('%s/%s_CRE%d.eps', chemin, prefix, q_method));

figure(3); plot(n_freq, if_hatmap(m_rg,t0+1)-lambda_0); grid;
xlabel('\lambda'); ylabel('IF estimation error [cycle/sample]')
title(sprintf('IF%d', if_method));
saveas(gcf, sprintf('%s/%s_IF%d.eps', chemin, prefix, if_method));

if_hatmap(m_0,t0+1)-lambda_0