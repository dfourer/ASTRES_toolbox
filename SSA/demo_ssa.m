% Automatic separation between a sinusoid and a chirp merged in a Gaussian white noise using SSA+CAH
% 
% Author: D. Fourer and J. Harmouche
% Date: Dec. 2016
% Contact: dominique@fourer.fr
%
% Refs: 
% [J. Harmouche, D. Fourer, P. Flandrin, F. Auger and P. Borgnat. One
% or Two Components ? The Singular Spectrum Analysis answers. Proc. SLRA'2015. Grenoble, France.June 2015]
% [J. Harmouche, D. Fourer, P. Flandrin, F. Auger and P. Borgnat. Une ou deux composantes:
% la réponse de l'analyse spectrale singulière. Proc. GRETSI'15. Lyon, France]

clear all
close all

N  = 300; %% signal length
L  = 40;  %% embedded dimension SSA parameter
nc = 2;   %% number of components

%% 1: Create signals
 time=(1:N)';
s=zeros(N,2);
% s1
lambda0=0.03; 
s(:,1) = cos(2*pi*lambda0*time);
s(:,1) = s(:,1) - mean(s(:,1));
% s2
lambda1     = 0.06;
deltalambda = (0.15-0.06);
s(:,2)  = 0.8 * cos(2*pi * (lambda1 * (time) + (deltalambda/(2*N)) * (time).^2));
s(:,2) = s(:,2) - mean(s(:,2));
% mix: x= s1 + s2
x = s(:,1)+s(:,2);              
% noise with arbitrary std value
noise = 0.05 * randn(N,1);
% mix + noise
xn = x+noise;


%% 2 Components extraction using SSA_auto
epsilon = 1e-3;    %% singular spectrum thresholding parameter (increase for more robustness to noise)
Y = ssa_decomp(xn, L, nc, epsilon);


%% 3 Display results
figure
plot(xn)
title(sprintf('observed mixture - SNR=%.2f dB', RQF(x, xn)))

figure
for i = 1:nc
 subplot(1,nc, i)
 plot(s(:,i))
 hold on
 plot(Y(:,i), 'r-.')
 legend('reference', 'reconstruction')
 title(sprintf('component %d - RQF=%.2f dB', i, RQF(s(:,i), Y(:,i))))
end

