%-------------------------------------------------------------------%
% 1D Decomposition algorithm allowing to extract the IMFs and trend %
%                                                                   %
% For theoretical aspects please refer to :                         %
%      http://arxiv.org/pdf/0911.1536                               %
%-------------------------------------------------------------------%


clear all
%close all
clc

addpath( genpath('include'), '-end' );
addpath( genpath('define'), '-end' );
addpath( genpath('EMDs'), '-end' );
addpath( genpath('package_emd'), '-end' );
addpath( genpath('emdos'), '-end' );
addpath( genpath('tftb-0.2'), '-end' );

%% Select an example
%type        = 'signal_2comp_ex1';
%type        = 'signal_2comp_ex2';
%type        = 'signal_2comp_ex3';
%type        = 'signal_2comp_ex4';
type            = 'signal_3comp_ex1';

%% Run algorithm
data            = create_signal(type);
algo            = define_parameters(type);
data.x0         = data.x;
results         = prox_emd(data.x, data, algo);


%% Print results
if data.K==2;
    figure(1); clf;
    subplot(2,2,1);plot(data.x0);title 'Original';
    subplot(2,2,2);plot(results(1,:));title 'Proposed EMD : IMF'
    subplot(2,2,4);plot(results(2,:));title 'Proposed EMD : trend';
elseif data.K==3
    figure(1); clf;
    subplot(3,2,1);plot(data.x0);title 'Original';
    subplot(3,2,2);plot(results(1,:));title 'Proposed EMD : IMF'
    subplot(3,2,4);plot(results(2,:));title 'Proposed EMD : IMF'   
    subplot(3,2,6);plot(results(3,:));title 'Proposed EMD : trend';   
end
