function sensitiv_io()
% Teste l'IO pour plusieurs methodes et k


% Signal
% Param
N = 1024*8;
T = 1;
t = linspace(0,T-T/N,N);
un = ones(size(t));


%% Signal s
a1 = 2*un+1.1*cos(4*pi*t);
phi1 = 2*pi*90;
f1=30*sin(3*pi*t);
fi1 = phi1+30*3*pi*cos(3*pi*t);

a2 = 1*un+0.4*sin(3*pi*t+0.2*un);
phi2 = 2*pi*49;
f2=21*sin(3*pi*t);
fi2 = phi2+21*3*pi*cos(3*pi*t);

a3 = 2*exp(-22*(t-T/2*un).^2);
phi3 = 2*pi*8*un;
a3 = 2*exp(-22*(t-T/2*un).^2);
phi3 = 2*pi*6*un;

s1 = a1.*sin(phi1*t+f1);
s2 = a2.*sin(phi2*t+f2);
s3 = a3.*sin(phi3.*t);
s = s1+s2+s3;%+sin(2*pi*180*t);
sred = s((N/4):(3*N/4));



% General options
opt.t = t;
opt.maxmodes = 2;
opt.stop = 'f';
opt.orderder = [2 2]; % Selection derivee
opt.liss = [0];



alpha = 0.02:0.05:0.3;
emd = zeros(size(alpha));
os4 = zeros(size(alpha));
os6 = zeros(size(alpha));
os8 = zeros(size(alpha));


opt.method = 'emd';
for i=1:length(alpha)
    opt.alpha = alpha(i);
    [imf,ort,nder] = emdos(s,opt);
    emd(i) = io(sred,imf(1:3,(N/4):(3*N/4)));
end



opt.method = 'os';
opt.sporder = 6;
for i=1:length(alpha)
    opt.alpha = alpha(i);
    [imf,ort,nder] = emdos(s,opt);disp(nder);
    os4(i) = io(sred,imf(1:3,(N/4):(3*N/4)));
end


opt.sporder = 8;
for i=1:length(alpha)
    opt.alpha = alpha(i);
    [imf,ort,nder] = emdos(s,opt);disp(nder);
    os6(i) = io(sred,imf(1:3,(N/4):(3*N/4)));
end


opt.sporder = 10;
for i=1:length(alpha)
    opt.alpha = alpha(i);
    [imf,ort,nder] = emdos(s,opt);disp(nder);
    os8(i) = io(sred,imf(1:3,(N/4):(3*N/4)));
end




% Affichage
figure();
plot(alpha,emd,'ko-',alpha,os4,'rs-',alpha,os6,'b*-',alpha,os8,'g+-');
legend('EMD','OS k=6','OS k=8','OS k=10');
xlabel('\alpha');ylabel('Orthogonality Index');
title('Orthogonality index for different spline orders');

