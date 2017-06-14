% test-script for VMD
% authors: Dominique Zosso and Konstantin Dragomiretskiy
% zosso@math.ucla.edu --- http://www.math.ucla.edu/~zosso
% Initial release 2013-12-12 (c) 2013
%
% When using this code, please do cite our paper:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing (in press)
% please check here for update reference: 
%          http://dx.doi.org/10.1109/TSP.2013.2288675

%--------------- Preparation
clear all;
close all;
clc;

% Time Domain 0 to T
T = 1000;
fs = 1/T;
t = (1:T)/T;
freqs = 2*pi*(t-0.5-1/T)/(fs);

% center frequencies of components
f_1 = 2;
f_2 = 24;
f_3 = 288;

% modes
v_1 = (cos(2*pi*f_1*t));
v_2 = 1/4*(cos(2*pi*f_2*t));
v_3 = 1/16*(cos(2*pi*f_3*t));

% for visualization purposes
fsub = {};
wsub = {};
fsub{1} = v_1;
fsub{2} = v_2;
fsub{3} = v_3;
wsub{1} = 2*pi*f_1;
wsub{2} = 2*pi*f_2;
wsub{3} = 2*pi*f_3;

% composite signal, including noise
f = v_1 + v_2 + v_3 + 0.1*randn(size(v_1));
f_hat = fftshift((fft(f)));

% some sample parameters for VMD
alpha = 2000;        % moderate bandwidth constraint
tau = 0;            % noise-tolerance (no strict fidelity enforcement)
K = 3;              % 3 modes
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-7;






%--------------- Run actual VMD code

[u, u_hat, omega] = VMD(f, alpha, tau, K, DC, init, tol);










%--------------- Visualization

% For convenience here: Order omegas increasingly and reindex u/u_hat
[~, sortIndex] = sort(omega(end,:));
omega = omega(:,sortIndex);
u_hat = u_hat(:,sortIndex);
u = u(sortIndex,:);
linestyles = {'b', 'g', 'm', 'c', 'c', 'r', 'k'};

figure('Name', 'Composite input signal' );
plot(t,f, 'k');
set(gca, 'XLim', [0 1]);

for sub = 1:length(fsub)
    figure('Name', ['Input signal component ' num2str(sub)] );
    plot(t,fsub{sub}, 'k');
    set(gca, 'XLim', [0 1]);
end

figure('Name', 'Input signal spectrum' );
loglog(freqs(T/2+1:end), abs(f_hat(T/2+1:end)), 'k');
set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
ylims = get(gca, 'YLim');
hold on;
for sub = 1:length(wsub)
    loglog([wsub{sub} wsub{sub}], ylims, 'k--');
end
set(gca, 'YLim', ylims);

figure('Name', 'Evolution of center frequencies omega');
for k=1:K
    semilogx(2*pi/fs*omega(:,k), 1:size(omega,1), linestyles{k});
    hold on;
end
set(gca, 'YLim', [1,size(omega,1)]);
set(gca, 'XLim', [2*pi,0.5*2*pi/fs], 'XGrid', 'on', 'XMinorGrid', 'on');

figure('Name', 'Spectral decomposition');
loglog(freqs(T/2+1:end), abs(f_hat(T/2+1:end)), 'k:');
set(gca, 'XLim', [1 T/2]*pi*2, 'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
hold on;
for k = 1:K
    loglog(freqs(T/2+1:end), abs(u_hat(T/2+1:end,k)), linestyles{k});
end
set(gca, 'YLim', ylims);


for k = 1:K
    figure('Name', ['Reconstructed mode ' num2str(K)]);
    plot(t,u(k,:), linestyles{k});   hold on;
    if ~isempty(fsub)
        plot(t, fsub{min(k,length(fsub))}, 'k:');
    end
    set(gca, 'XLim', [0 1]);
end

