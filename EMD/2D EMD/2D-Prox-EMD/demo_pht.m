%-------------------------------------------------------------------%
% 2D Prony-Huang Transform : spectral analysis of an image          %
% combining EMD decomposition and local Prony spectral estimation   %
%                                                                   %
% For theoretical aspects please refer to :                         %
%      http://arxiv.org/pdf/1404.7680v1.pdf                         %
%-------------------------------------------------------------------%

clear all
close all
clc

addpath( genpath('include'), '-end');
addpath( genpath('colormap_utilities'), '-end');

%select an example
example = 'example1';
%example = 'example2';

%generate data

data = create_signal(example);
data.x0 = data.x;

%Chose type of 2D EMD

%type = 'G2D';          %Genuine 2D
%type = 'G2Ddir';       %Directional genuine 2D
%type = 'P2DLC';        %Pseudo 2D on lines and columns
type = 'P2DLCD';        %Pseudo 2D on lines, columns, diagonals and anti-diagonals

%define parameters

lambda = define_parameters_emd_2D(example,type); % parameters of EMD algorithm
sizepatch = define_sizepatch(example);           % size of patches for Prony spectral analysis

%run EMD algorithm

resultsemd = run_prox_emd_2D(data,type,lambda);

%run Prony algorithm

resultsprony = run_prony_2D(resultsemd.imf,sizepatch,lambda.K);

%plot results
figure(1);
subplot(1,lambda.K + 1,1);imagesc(resultsemd.trend);axis image off; colormap(gray);
title 'Trend'
for indi=1:lambda.K
    subplot(1,lambda.K+1,indi+1);imagesc(resultsemd.imf{indi});axis image off; colormap(gray);
    tit=['IMF ',num2str(indi)];
    title(tit)
end

figure(2);
for indi=1:lambda.K
    subplot(1,lambda.K,indi);semilogy(resultsemd.crit{indi}); hold on;
    tit=['crit ',num2str(indi)];
    title(tit)
end

for indi=1:lambda.K
    figure(2+indi);
    subplot(221);imagesc(resultsemd.imf{indi});axis image off;colormap(gray);
    title 'IMF'
    freezeColors;
    cbfreeze(colorbar);
    subplot(222);imagesc(resultsprony.amplitude{indi},[0 max(max(resultsprony.amplitude{indi}))]);axis image off;colormap(gray);
    title 'amplitude'
    freezeColors;
    cbfreeze(colorbar);
    subplot(223);imagesc(resultsprony.freq{indi},[0 1]);axis image off;colormap(hsv);colorbar;
    title 'frequency'
    freezeColors;
    cbfreeze(colorbar);
    subplot(224);imagesc(resultsprony.orientation{indi},[-pi/2 pi/2]);axis image off;colormap(hsv);colorbar;
    title 'orientation'
    freezeColors;
    cbfreeze(colorbar);
end