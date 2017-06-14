function scfig2()
% scfig 2 : plots figure 2 of the paper "On the mode synthesis in the
% synchrosqueezing method"


[s s1 s2 s3] = gentests();
N = 1024;
dt = 1/N;
t = dt*(0:N-1);
nv = 32;
mywav = 'cmor8-1';
nr = 3;
lambda = 100;
clwin = 2;
gamma = 0.01;
reg = 0;
        
%s = sin(2*pi*50*t);

% Wavelet
[Wx Tx fs as] = sqt(s,dt,gamma,mywav,nv);

% Ridge extraction
[Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin);
% Computating a posteriori the wavelet mask SG
Delta = floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps);
centpsi = get_centpsi(mywav);

% Reconstruction, old method
imf1 = synth_tx(Tx,dt,mywav,nv,Cs,clwin);
SG1 = compute_mask(Wx,fs,Cs,clwin,dt,gamma);

% Recondstruction new method
SG2 = compute_maskWx(Wx,Cs,Delta,centpsi,fs);
SG2 = flipud(SG2);
imf2 = synth_mult(Wx,dt,mywav,nv,SG2,0,as);

% Displaying TX, masks and Waveket transform
figure();ha = gca();sqplot(flipud(exp(SG1-1)),dt,fs,ha);
colormap('gray');set(gcf,'colormap',flipud(get(gcf,'colormap')));
figure();ha = gca();sqplot(flipud(exp(SG2-1)),dt,fs,ha);
colormap('gray');set(gcf,'colormap',flipud(get(gcf,'colormap')));
figure();ha = gca();sqplot(flipud(Wx),dt,fs,ha);
colormap('gray');set(gcf,'colormap',flipud(get(gcf,'colormap')));
figure();ha = gca();sqplot(Tx,dt,fs,ha);
colormap('gray');set(gcf,'colormap',flipud(get(gcf,'colormap')));

% Display
figure();
plot(t,imf1(3,:),'b-','LineWidth',2);
hold on;
plot(t,s2,'r--','LineWidth',2);
xlabel('t');
legend('Old reconstruction','Real mode');%,'Meignen');
set(gca,'xlim',[0.15 0.35]);
%return;

% Display
figure();
plot(t,imf2(3,:),'b-','LineWidth',2);
hold on;
plot(t,s2,'r--','LineWidth',2);
xlabel('t');
legend('New reconstruction','Real mode');%,'Meignen');
set(gca,'xlim',[0.15 0.35]);

disp('');
