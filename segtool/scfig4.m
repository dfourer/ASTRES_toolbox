function scfig4()
% scfig 4 : plots figure 4 of the paper "On the mode synthesis in the
% synchrosqueezing method"

[s s1 s2 s3] = gentests();
N = 1024;
dt = 1/N;
t = dt*(0:N-1);
nv = 32;
mywav = 'cmor8-1';
nr = 3;
lambda = 100;
clwin = 3;
gamma = 0.01;
reg = 0;
        


x = awgn(s,2);

%% NO NOISE : getting the right curve c
[Wx Tx fs as] = sqt(s,dt,gamma,mywav,nv);
% Ridge extraction
[Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin);
% Computating a posteriori the wavelet mask SG
Delta = floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps);
centpsi = get_centpsi(mywav);
% Perturbation
dk = floor(N/30);
for k=floor(N/10):floor(N/10):floor(N-N/6);
    Cs(3,k:k+dk) = Cs(3,k);
end
SG = compute_maskWx(Wx,Cs,Delta,centpsi,fs);
SG = flipud(SG);

% Displaying mask
figure();ha = gca();sqplot(flipud(exp(SG-1)),dt,fs,ha);


%% NOISY signal
[Wx Tx fs as] = sqt(x,dt,gamma,mywav,nv);

figure();ha = gca();sqplot(flipud(Wx),dt,fs,ha);

SG(SG==1)=0;SG(SG==2)=0;SG(SG==3)=1;
% Reconstruction, without regularization
imf1 = synth_mult(Wx,dt,mywav,nv,SG,0,as);
% Display
figure();
plot(t,imf1,'b-','LineWidth',2);
hold on;
plot(t,s2,'r--','LineWidth',2);
xlabel('t');
legend('Old reconstruction','Real mode');%,'Meignen');
set(gca,'xlim',[0.38 0.48]);
%return;

% Reconstruction with regularization
imf3 = synth_mult(Wx,dt,mywav,nv,SG,1,as); % RKHS regularization
% Display
figure();
plot(t,imf3,'b-','LineWidth',2);
hold on;
plot(t,s2,'r--','LineWidth',2);
xlabel('t');
legend('New reconstruction','Real mode');%,'Meignen');
set(gca,'xlim',[0.38 0.48]);


return;

%% Old method with coefficients set X ///  NOT USED 
imf0 = synth_tx(Tx,dt,mywav,nv,Cs,clwin);
SG0 = compute_mask(Wx,fs,Cs,clwin,dt,gamma);
W = zeros(size(Wx));W(SG0==3)=Wx(SG0==3);
figure();ha = gca();sqplot(flipud(W),dt,fs,ha);
colormap('gray');set(gcf,'colormap',flipud(get(gcf,'colormap')));
% Display
figure();
plot(t,imf0(3,:),'b-','LineWidth',2);
hold on;
plot(t,s2,'r--','LineWidth',2);
xlabel('t');
legend('old reconstruction set X','Real mode');%,'Meignen');
set(gca,'xlim',[0.38 0.48],'ylim',[-1.5 2]);


