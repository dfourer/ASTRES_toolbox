clear all
%close all
clc

addpath( genpath('include'), '-end' );
addpath( genpath('define'), '-end' );
addpath( genpath('EMDs'), '-end' );
addpath( genpath('package_emd'), '-end' );
addpath( genpath('emdos'), '-end' );
addpath( genpath('tftb-0.2'), '-end' );

%% Create signal
N          = 500;   t = kron(ones(1,N),[1:-0.1:-1]); 
f          = [0.2, 0.15,0.01] ;
df         = [0.05, 0.05, 0.05];
a          = [2, 1, 0.5];
x1         = a(1).*real(fmsin(N, f(1)-df(1), f(1)+df(1), N));
x2         = amgauss(N,N/2,4*N/3).*real(fmlin(N, f(2)-df(2), f(2)+df(2), N));
x3         = a(3).*real(fmconst(N, f(3))).*([zeros(1,125), ones(1,200),zeros(1,100),ones(1,75)])';% ;triangular_signal(tmp.N,20)';%triangular_signal(tmp.N,55);
x          = x1' +x2' + x3';


%% Run Prox-EMD
%results     = prox_emd(x);

%% Run Prox-EMD with options
data.K           = 3;
algo.nbiter      = [3000, 10000];
algo.epsa        = [40,0.05];
algo.epsd        = [0.05,0.001];
algo.dataFid{1}  = 'l2';
algo.dataFid{2}  = 'l2';
algo.regTre{1}   = 'l1';
algo.regTre{2}   = 'l2';
algo.filt{1}     = 'laplacian';    
algo.filt{2}     = 'laplacian'; 
algo.regFlu{1}   = 'l1';    
algo.regFlu{2}   = 'l1'; 
algo.ext{1}      = 'd1'; 
algo.ext{2}      = 'd1'; 
results          = prox_emd(x, data, algo);


%% Print results
figure(1);clf; 
subplot(3,3,1);plot(x);title 'Original';axis([0 length(x) 1.1*min(x) 1.1*max(x)]);
subplot(3,3,2);plot(x1);title 'Original: IMF 1'; axis([0 length(x) 1.1*min(x1) 1.1*max(x1)]);
subplot(3,3,5);plot(x2);title 'Original:  IMF 2';  axis([0 length(x) 1.1*min(x2) 1.1*max(x2)]); 
subplot(3,3,8);plot(x3);title 'Original:  Trend'; axis([0 length(x) 1.1*min(x3) 1.1*max(x3)]); 
subplot(3,3,3);plot(results(1,:));title 'Proposed EMD: IMF 1'; axis([0 length(x) 1.1*min(x1) 1.1*max(x1)]);
subplot(3,3,6);plot(results(2,:));title 'Proposed EMD:  IMF 2' ; axis([0 length(x) 1.1*min(x2) 1.1*max(x2)]); 
subplot(3,3,9);plot(results(3,:));title 'Proposed EMD:  Trend'; axis([0 length(x) 1.1*min(x3) 1.1*max(x3)]);  