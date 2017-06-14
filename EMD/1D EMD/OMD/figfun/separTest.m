function separTest()
%% separTest : test of separation power, comparison with the original EMD

%% Paramètres
N = 1048*8;
T = 32;
t = linspace(0,T,N);
s1 = sin(2*pi*t);
intval = N/4:3*N/4; % remove border effects

% General options
opt.t = t;
opt.MAXMODES = 1;
opt.alpha= 0.01;
opt.sporder=10;
opt.orderder = -1;

f = 0.01:0.01:0.99;
%f = 0.01:0.1:0.99;
P = length(f);

a = 10.^(-2:0.025:2);
%a = 10.^(-2:0.5:-0.5);
M = length(a);

ordOS = zeros(P,M); % OS
ordEMD = zeros(P,M); % EMD
ordEMDNI = zeros(P,M); % EMDNI


for i=1:P
    for j=1:M
        s2 = a(j)*sin(2*pi*f(i)*t);
        %s2 = 0.2*sin(2*pi*0.3*t);
        s = s1+s2;
        
        % OS
        opt.method='os';
        imf = emdos(s,opt);
        imf1=imf(1,:);
        ordOS(i,j) = norm(imf1(intval)-s1(intval))/norm(s2(intval));
        %figure();plot_imf(imf,t,'');title('os');        disp(ordOS(i,j));
        
        % EMDNI
        opt.method='emdni';
        imf = emdos(s,opt);
        imf1=imf(1,:);
        ordEMDNI(i,j) = norm(imf1(intval)-s1(intval))/norm(s2(intval));
        
        % EMD
        opt.method='emd';
        imf = emdos(s,opt);
        imf1=imf(1,:);
        ordEMD(i,j) = norm(imf1(intval)-s1(intval))/norm(s2(intval));
        %figure();plot_imf(imf,t,'');disp(ordOS(i,j));title('emd');
        
     end
end

xx = -2:0.05:2;

figure();
imagesc(log(a)/log(10),f,ordOS);
set(gca,'YDir','normal')
colormap('gray');
xlabel('a (log_{10})');
ylabel('f');
zlabel('erreur normalisee (log)');
title('Separation power, OS, automatic selection of derivation order');
h=findall(gcf,'type','line');delete(h);h=findall(0,'type','legend');delete(h);set(gca,'CLim',[0 1.1]);hold on;plot(xx,10.^(-1/1*xx),'color',0.9*ones(1,3),'LineWidth',2);ylim([0 1]);xlim([-2 2]);plot(xx,10.^(-1/2*xx),'color',0.7*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/3*xx),'color',0.4*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/4*xx),'color',0.1*ones(1,3),'LineWidth',2);legend('af=1','af^2=1','af^3=1','af^4=1');
colorbar;

figure();
imagesc(log(a)/log(10),f,ordEMD);
set(gca,'YDir','normal')
colormap('gray');
xlabel('a (log_{10})');
ylabel('f');
zlabel('erreur normalisee (log)');
title('Separation power, Standard EMD with stopping criteria');
h=findall(gcf,'type','line');delete(h);h=findall(0,'type','legend');delete(h);set(gca,'CLim',[0 1.1]);hold on;plot(xx,10.^(-1/1*xx),'color',0.9*ones(1,3),'LineWidth',2);ylim([0 1]);xlim([-2 2]);plot(xx,10.^(-1/2*xx),'color',0.7*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/3*xx),'color',0.4*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/4*xx),'color',0.1*ones(1,3),'LineWidth',2);legend('af=1','af^2=1','af^3=1','af^4=1');
colorbar;

figure();
imagesc(log(a)/log(10),f,ordEMDNI);
set(gca,'YDir','normal')
colormap('gray');
xlabel('a (log_{10})');
ylabel('f');
zlabel('erreur normalisee (log)');
title('Separation power, EMD-NI');
h=findall(gcf,'type','line');delete(h);h=findall(0,'type','legend');delete(h);set(gca,'CLim',[0 1.1]);hold on;plot(xx,10.^(-1/1*xx),'color',0.9*ones(1,3),'LineWidth',2);ylim([0 1]);xlim([-2 2]);plot(xx,10.^(-1/2*xx),'color',0.7*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/3*xx),'color',0.4*ones(1,3),'LineWidth',2);plot(xx,10.^(-1/4*xx),'color',0.1*ones(1,3),'LineWidth',2);legend('af=1','af^2=1','af^3=1','af^4=1');
colorbar;


end