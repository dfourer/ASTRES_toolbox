function sensitiv_spo()
% Test function, creates figure 4 of [1].


N = 4*1024;
[s s1 s2 s3] = gen_tests();


figure(1);hold on;
figure(2);hold on;

p_alpha = [0.02:0.003:0.06 0.06:0.01:0.2];
index = N/4+1:3*N/4; % remove boundary effects


%% EMD

% General options
opt.maxmodes = 1;
opt.sporder = 8;
opt.alpha= 0.2;
opt.stop = 'f';
opt.orderder = [0 0 0]; % Selection derivee
opt.method = 'emdni';
opt.liss = [0];



err = zeros(size(p_alpha));
errp = zeros(size(p_alpha));
oldimf=zeros(1,N);
for j=1:length(p_alpha)
    opt.alpha = p_alpha(j);
    [imf,ort,nder] = emdos(s,opt);
    errp(j) = norm(oldimf(1,index)-imf(1,index));
    oldimf = imf;
    err(j) = norm(s1(index)-imf(1,index));
end

figure(1);
plot(p_alpha,log(err*(1/norm(s))),'ko--');
figure(2);
plot(p_alpha(2:end),1/sqrt(N/2)*(errp(2:end)./diff(p_alpha)),'ko--');



%% EMDOS

% General options
opt.maxmodes = 1;
opt.sporder = 8;
opt.alpha= 0.06;
opt.stop = 'f';
opt.orderder = [2 0 0]; % Selection derivee
opt.method = 'os';
opt.liss = [0];



err6 = zeros(size(p_alpha));
errp6 = zeros(size(p_alpha));
err8 = zeros(size(p_alpha));
errp8 = zeros(size(p_alpha));
err10 = zeros(size(p_alpha));
errp10 = zeros(size(p_alpha));
oldimf6=zeros(1,N);
oldimf8=zeros(1,N);
oldimf10=zeros(1,N);
for j=1:length(p_alpha)
    opt.alpha = p_alpha(j);
    % k=6
    opt.sporder=6;
    [imf,ort,nder] = emdos(s,opt);
    errp6(j) = norm(oldimf6(1,index)-imf(1,index));
    oldimf6 = imf;
    err6(j) = norm(s1(index)-imf(1,index));
    % k=8
    opt.sporder=8;
    [imf,ort,nder] = emdos(s,opt);
    errp8(j) = norm(oldimf8(1,index)-imf(1,index));
    oldimf8 = imf;
    err8(j) = norm(s1(index)-imf(1,index));
    % k=6
    opt.sporder=10;
    [imf,ort,nder] = emdos(s,opt);
    errp10(j) = norm(oldimf10(1,index)-imf(1,index));
    oldimf10 = imf;
    err10(j) = norm(s1(index)-imf(1,index));
end

figure(1);
plot(p_alpha,log(err6/norm(s)),'rs-');
plot(p_alpha,log(err8/norm(s)),'b*-');
plot(p_alpha,log(err10/norm(s)),'g+-');
figure(2);
plot(p_alpha(2:end),1/sqrt(N/2)*(errp6(2:end)./diff(p_alpha)),'rs-');
plot(p_alpha(2:end),1/sqrt(N/2)*(errp8(2:end)./diff(p_alpha)),'b*-');
plot(p_alpha(2:end),1/sqrt(N/2)*(errp10(2:end)./diff(p_alpha)),'g+-');




%% MISE EN FORME
figure(1);
legend('EMD','OS k=6','OS k=8','OS k=10');
xlabel('\alpha');
set(findall(0,'type','text'),'FontSize',15);
set(gca,'FontSize',15);
set(findall(0,'type','line'),'linewidth',2);

figure(2);
legend('EMD','OS k=6','OS k=8','OS k=10');
xlabel('\alpha');
set(findall(0,'type','text'),'FontSize',15);
set(gca,'FontSize',15);
set(findall(0,'type','line'),'linewidth',2);


end