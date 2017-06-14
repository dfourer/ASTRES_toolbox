function [SG Nf imf] = ts_segment(s,Wx,as,dt,mywav,nv,nmod,rrec,gamm,reg)
% Function ts_segment : segments the time-scale plane according to paper of
% Meignen, Oberlin, and McLaughlin.


% Setting the extremal values of test set Gamma
if strcmp(gamm,'auto')
    sigma = estim_sigma(s);
    Gammin = 1/sigma/sqrt(2)/2;
else
    Gammin = str2num(gamm);
end

% Other prameters
Ndis = 100;
[na N] = size(Wx);

% Estimates the number of modes
if strcmp(nmod,'auto')
    Nf = estim_nf(Wx,as,Gammin,Ndis);
else
    Nf = str2num(nmod);
end


Delta = floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps);
% Detecting intervals
intj = zeros(2,Nf);
for b=1:N % Time loop
    vec = abs(Wx(:,b));
    [flag gam] = get_gam(vec,Nf,Ndis,Gammin);
    if ~flag
        disp('warning! not the right number of modes!!');
        %errordlg('Did not manage to recover the right number of modes');
        %return
    else
        cpt = 1;
        for a=1:na-1
            if vec(a)<gam && vec(a+1)>=gam
                intj(1,cpt) = a;
            end
            if vec(a)>=gam && vec(a+1)<gam
                intj(2,cpt) = a;
                cpt = cpt+1;
            end
        end
        % Center and size of support
        SG(:,b) = mplex(floor(0.5*(intj(1,:)+intj(2,:))),Nf,Delta,as);
    end
end

% Computing IMFs
imf = synth_mult(Wx,dt,mywav,nv,SG,reg);

return;


%% RKHS experiment
% a = floor(0.4*na);
% b = floor(0.25*N);
% 
% tic
% %Kab = rkhs(mywav,dt,as,N,a,b); % Just Reproducing kernel
% Wtilde = projrkhs(Wx,mywav,dt,as,N,nv);
% toc
% 
% 
% %figure();imagesc(abs(Kab));figure();imagesc(angle(Kab));
% %figure();imagesc(log(1+abs(Wtilde)));figure();imagesc(log(1+abs(Wx)));
% 
