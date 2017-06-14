function scfig3()
% scfig 3 : plots figure 3 of the paper "On the mode synthesis in the
% synchrosqueezing method"

[s s1 s2 s3] = gentests();
P = 120;
N = 1024;
dt = 1/N;

nv = 32;
mywav = 'cmor8-1';
nr = 3;
lambda = 100;
clwin = 3;
gamma = 0.01;
reg = 0;
        
snr = linspace(-3,10,P);
err1 = zeros(size(snr));
err2 = zeros(size(snr));
err3 = zeros(size(snr));

for k=1:P
    % Adding noise
    x = awgn(s,snr(k));
    
    [Wx Tx fs as] = sqt(x,dt,gamma,mywav,nv);
    % Ridge extraction
    [Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin);
    
    % Computating a posteriori the wavelet mask SG
    Delta = floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps);
    centpsi = get_centpsi(mywav);
    SG = compute_maskWx(Wx,Cs,Delta,centpsi,fs);
    SG = flipud(SG);
    imf1 = synth_mult(Wx,dt,mywav,nv,SG,reg,as);
    imf2 = synth_tx(Tx,dt,mywav,nv,Cs,clwin);
    
    
    err1(k) = norm(imf1(1,:)-s3);
    err2(k) = norm(imf2(1,:)-s3);
end

err1 = 10*log(err1/norm(s2));
err2 = 10*log(err2/norm(s2));

% Display
figure();
plot(snr,err2,'b-','LineWidth',2);
hold on;
plot(snr,err1,'r--','LineWidth',2);
plot(snr,err3,'k--','LineWidth',2);
xlabel('SNR (dB)');
ylabel('MSE (dB)');
legend('Old reconstruction','New reconstruction');%,'Meignen');