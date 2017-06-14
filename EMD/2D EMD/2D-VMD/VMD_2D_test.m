% Test script for 2D-VMD
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% {konstantin,zosso}@math.ucla.edu
% http://www.math.ucla.edu/~{konstantin,zosso}
% Initial release 2014-03-17 (c) 2014
%
% When using this code, please do cite our papers:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing, 62(3):531-544, 2014. DOI:10.1109/TSP.2013.2288675
%
% K. Dragomiretskiy, D. Zosso, Two-Dimensional Variational Mode
% Decomposition, IEEE Int. Conf. Image Proc. (submitted). Preprint
% available here: ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf
%

%% preparations

close all;
clc;
clear all;


% Sample data
texture = load('texture.mat');
f = texture.f;

% parameters:
alpha = 5000;       % bandwidth constraint
tau = 0.25;         % Lagrangian multipliers dual ascent time step
K = 5;              % number of modes
DC = 1;             % includes DC part (first mode at DC)

init = 1;           % initialize omegas randomly, may need multiple runs!

tol = K*10^-6;      % tolerance (for convergence)

%% run actual 2D VMD code

[u, u_hat, omega] = VMD_2D(f, alpha, tau, K, DC, init, tol);


%% Visualization
figure('Name', 'Input image');
imagesc(f);
colormap gray;
axis equal;
axis off;

figure('Name', 'Input spectrum');
imagesc(abs(fftshift(fft2(f))));
colormap gray;
axis equal;
axis off;
hold on;
colors = 'rbycmg';
for k = 1:size(omega,3)
    plot( size(f,2)*(0.5+omega(:,1,k)), size(f,1)*(0.5+omega(:,2,k)), colors(k) );
end

for k=1:size(u,3)
    figure('Name', ['Mode #' num2str(k)]);
    imagesc(u(:,:,k));
    colormap gray;
    axis equal;
    axis off;
end

figure('Name', 'Reconstructed composite');
imagesc(sum(u,3));
colormap gray;
axis equal;
axis off;

