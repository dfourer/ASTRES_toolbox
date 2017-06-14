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
example = 'example2';

%generate data

data = create_signal(example);
data.x0 = data.x;

%Chose type of 2D EMD

type = 'G2D';          %Genuine 2D
%type = 'G2Ddir';       %Directional genuine 2D
%type = 'P2DLC';        %Pseudo 2D on lines and columns
%type = 'P2DLCD';        %Pseudo 2D on lines, columns, diagonals and anti-diagonals

wave3 = generate_wave(10,40,-pi/4,256);


%define parameters
i=1;
j=1;

for lambda1 = 0.01
    for lambda2 = [15]

lambda = define_parameters_emd_2D_lambda(example,type,lambda1,lambda2);

%run algorithm

results = run_prox_emd_2D(data,type,lambda);
error(i,j) = sqrt(norm(results.imf{1} - wave3,'fro'));
j=j+1;
    end
    i=1+i;
    j=1;
end

%plot results

% figure(1);hold on;
% map=colormap(hsv(3));
% for i = 1:3
%     lambda2=[1 5 10 15 20];
%    plot(lambda2,error(i,:),'Color',map(i,:));
% end
% hold off;

figure(1);
%map=colormap(hsv(3));
    lambda2=[1 5 10 15 20];
   plot(lambda2,error(1,:));


% %plot results
% figure(1);
% subplot(1,lambda.K + 1,1);imagesc(results.trend);axis image off; colormap(gray);
% title 'Trend'
% for indi=1:lambda.K
%     subplot(1,lambda.K+1,indi+1);imagesc(results.imf{indi});axis image off; colormap(gray);
%     tit=['IMF ',num2str(indi)];
%     title(tit)
% end
% 
% figure(2);
% for indi=1:lambda.K
%     subplot(1,lambda.K,indi);semilogy(results.crit{indi}); hold on;
%     tit=['crit ',num2str(indi)];
%     title(tit)
% end
