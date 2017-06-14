clear all
close all

addpath('./ST')

N = 500;
mu = 1e-4;
x = 0.01 * randn(1, N);
t = 1:N;

Nh = 301; %81; %127;% short-time window length
m_max = N/2;

w = tftb_window(Nh,'Kaiser');
%[sp,rs] = tfrrsp(x,1:length(x),Nf,w);
%[tfr,a,b] = tfrlsp(x(:),t,N,w);

seuil = 1e-2;

[tfr,rtfr,hat] = tfrlmrgab(x(:), t, N, mu, Nh);

R_n = repmat(1:N, N,1); %% time
R_m = repmat((1:N).', 1,N);  %% frequency
r = 1j * R_n + R_m;
%r = zeros(size(tfr));

[I,J] = find(abs(hat-r) < seuil);

nb_z = length(I);

figure(979)
imagesc(abs(tfr))
colormap gray;
colormap(flipud(colormap));
hold on

plot(I, J, 'rx');

%%%%%%%%%%%%%%%%%%%

x2 = x + real(fmconst(N, 0.2)).';
[tfr2,rtfr2,hat2] = tfrlmrgab(x2(:), t, N, mu, Nh);

[I2,J2] = find(abs(hat2-r) < seuil);

figure(980)
imagesc(abs(tfr2))
colormap gray;
colormap(flipud(colormap));
hold on
plot(I2, J2, 'rx');
 
figure
plot(I, J, 'rx');
hold on
plot(I2, J2, 'g.');




