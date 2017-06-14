function scfig1()
% scfig 1 : plots figure 1 of the paper "On the mode synthesis in the
% synchrosqueezing method"


N = 1024;
dt = 1/1024;
t = dt*(0:N-1);

s = cos(2*pi*15*t)+cos(2*pi*25*t);
mywav = 'cmor8-1';
nv = 32;

% Wavelet transform
[Wx Tx fs as] = sqt(s,dt,0.01,mywav,nv,0);
%h=figure();ha = gca(h);sqplot(Wx,dt,as,ha);return;

na = size(Wx,1);

% Selecting a ridge
W = zeros(size(Wx));
delt = 6;

% Piecewise constant 
idx = [floor(0.6*na)*ones(1,N/4) floor(0.7*na)*ones(1,N/4) floor(0.6*na)*ones(1,N/4) floor(0.7*na)*ones(1,N/4)];
% Modulated, smooth
idx = floor(na*(0.52+0.05*cos(4*pi*t)));

for k=1:N
    W(idx(k)+(-delt:delt),k) = Wx(idx(k)+(-delt:delt),k);
end
% Display
h=figure();ha = gca(h);
sqplot(flipud(W),dt,as,ha);%return;
xlabel('t');ylabel('1/a');
set(gca,'xlim',[0.1 0.9],'ylim',[100 400]);
colormap('gray');

% synthesis
imf = synth_sq(W,dt,mywav,nv);
figure();plot(t,imf(1,:));

% Wavelet transform and display
Wnew = mycwt(imf,as,mywav,dt);
h=figure();ha = gca(h);
sqplot(flipud(Wnew),dt,as,ha);
set(gca,'xlim',[0.1 0.9],'ylim',[100 400]);
colormap('gray');