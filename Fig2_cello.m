clear all
close all

chemin = './figs';
if ~exist(chemin, 'dir')
  mkdir(chemin);
end

addpath('./EMD/1D EMD/Classical EMD')
addpath('./ST')
addpath('./SSA')
addpath('./SSA')
addpath('./segtool')

signal_name = 'Cello';
[x, Fs] = wavread('cello_short.wav');
x = x(1:1600);
%x = x(1:1000);

%% tfr representation
Nh = 301; %81; %127;% short-time window length
Nf = 1024; %256;% # of frequency bins
m_max = Nf/2;

w = tftb_window(Nh,'Kaiser');
%[sp,rs] = tfrrsp(x,1:length(x),Nf,w);
Fx     = tfrstft(x,1:length(x),Nf,w);
sp = abs(Fx).^2;
f = Fs*((1:m_max)-1)/Nf;
t = ((1:length(x))-1)/Fs;
figure(21)
% subplot(121)
% imagesc(t, f, rs(1:m_max,:).^0.15); %1:128  % 60/2048 * Fe1 = 480Hz
% set(gca,'Ydir','normal');
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% title('reassigned spectrogram')
%subplot(122)

imagesc(t, f, sp(1:m_max,:).^0.15);
set(gca,'Ydir','normal');
colormap gray;
colormap(flipud(colormap));
xlabel('time [s]')
ylabel('frequency [Hz]')
title(sprintf('%s spectrogram', signal_name));
saveas(gcf, sprintf('%s/cello_spectrogram.eps', chemin), 'epsc');
pause(1)


%% band-pass filtering
% mask = zeros(size(Fx));
% m_range = 60:320;
% mask([m_range Nf-m_range],:) = 1;
% Fx_filtered = Fx .* mask;
% figure
% imagesc(t, f, abs(Fx_filtered(1:m_max,:)).^0.15); %1:128  % 60/2048 * Fe1 = 480Hz
% set(gca,'Ydir','normal');
% xlabel('time [s]')
% ylabel('frequency [Hz]')
% title('filtered signal spectrogram')
% pause(1)
% xref = x;
% [x] = real(tfristft(Fx_filtered, 1:length(x),w));


%% 2 SSA
nc     = 4;
L      = 150;  %% frequency resolution 1/L * Fs
epsilon = 6e-2;

%% test classic ssa
Y2 = ssa_decomp(x, L, nc, epsilon);
[ Y2 ] = sort_components( Y2 );
%plot_comp(Y2, [], 0, [-0.5 0.5]);

figure(231)
for i = 1:nc
 subplot(nc,1,i)
 [sp] = tfrsp(Y2(:,i), 1:length(Y2(:,i)), Nf, w);
 imagesc(t, f, sp(1:m_max,:).^0.15) %1:128  % 60/2048 * Fe1 = 480Hz
 set(gca,'Ydir','normal');
 colormap gray;colormap(flipud(colormap));
 if i == nc,  xlabel('time [s]'); end
 ylabel('frequency [Hz]')
 if i == 1, title('SSA extracted components'); end
end
saveas(gcf, sprintf('%s/cello_SSA.eps', chemin), 'epsc');
pause(1)


%% test EMD
Y_emd = emd(x, 'MAXMODES', nc);
Y_emd = Y_emd.';
[ Y_emd ] = sort_components( Y_emd);
%plot_comp(Y_emd(:,1:4), [], 0, [-0.5 0.5]);

figure(232)
for i = 1:nc
 subplot(nc,1,i)
 [sp] = tfrsp(Y_emd(:,i), 1:length(Y2(:,i)), Nf, w);
 imagesc(t, f, sp(1:m_max,:).^0.15) %1:128  % 60/2048 * Fe1 = 480Hz
 set(gca,'Ydir','normal');
 colormap gray;colormap(flipud(colormap));
 if i == nc,  xlabel('time [s]'); end
 ylabel('frequency [Hz]')
 if i == 1, title('EMD extracted components'); end
end
saveas(gcf, sprintf('%s/cello_EMD.eps', chemin), 'epsc');
pause(1)

%% test synchrosqueezing
M = 1024;
L = 20;
gamma_K = 1e-4;
[tfr, stfr] = my_tfrsgab(x, M, L, gamma_K);     %% (faster) classical  synchrosqueezed STFT

%% ridge extraction
Mh = round(M/2);
m_range = m_axis(M)/M;
lambda = 0.01; %1e-1;    %% ridge detection parameter 1
clwin  = 9;              %% ridge detection parameter 2
K      = 3;              %% vicinity of the ridge

[ mask ] = ridge_detect_brvmask(stfr(1:Mh,:), m_range(1:Mh), nc, lambda, clwin, K);

for i = 1:nc
  mask2(:,:,i) = [mask(:,:,i);mask(end:-1:1,:,i)];
end
mask = mask2;

debug_trace = 0;
[ Y_sst ] = sst_comp_ext(stfr, mask, L, debug_trace);
Y_sst = Y_sst.';
[ Y_sst ] = sort_components( Y_sst);
%plot_comp(real(Y_sst), [], 0, [-0.5 0.5]);

figure(234)
for i = 1:nc
 subplot(nc,1,i)
 [sp] = tfrsp(Y_sst(:,i), 1:length(Y2(:,i)), Nf, w);
 imagesc(t, f, sp(1:m_max,:).^0.15) %1:128  % 60/2048 * Fe1 = 480Hz
 set(gca,'Ydir','normal');
 colormap gray;colormap(flipud(colormap));
 if i == nc,  xlabel('time [s]'); end
 ylabel('frequency [Hz]')
 if i == 1, title('synchrosqueezing extracted components'); end
end
saveas(gcf, sprintf('%s/cello_synchrosqueezing.eps', chemin), 'epsc');
pause(1)
eps2pdf(chemin);