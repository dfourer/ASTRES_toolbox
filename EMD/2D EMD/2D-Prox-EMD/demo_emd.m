%-------------------------------------------------------------------%
% 2D variational EMD : decomposition of an image into trend and imf                                             %
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
%example = 'example1';
%example = 'example2';
example = 'example3';


%generate data

data = create_signal(example);
data.x0 = data.x;

%Chose type of 2D EMD

type = 'G2D';          %Genuine 2D
%type = 'G2Ddir';       %Directional genuine 2D
%type = 'P2DLC';        %Pseudo 2D on lines and columns
%type = 'P2DLCD';        %Pseudo 2D on lines, columns, diagonals and anti-diagonals

%define parameters

lambda = define_parameters_emd_2D(example,type);

%run algorithm

results = run_prox_emd_2D(data,type,lambda);

%plot results
figure(1);
subplot(1,lambda.K + 1,1);imagesc(results.trend);axis image off; colormap(gray);
title 'Trend'
for indi=1:lambda.K
    subplot(1,lambda.K+1,indi+1);imagesc(results.imf{indi});axis image off; colormap(gray);
    tit=['IMF ',num2str(indi)];
    title(tit)
end

figure(2);
for indi=1:lambda.K
    subplot(1,lambda.K,indi);semilogy(results.crit{indi}); hold on;
    tit=['crit ',num2str(indi)];
    title(tit)
end
