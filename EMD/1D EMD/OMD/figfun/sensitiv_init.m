function sensitiv_init()
% creates figure 5 of [1] : evaluates the sensitivity of the method to the
% initial guess \hat m.
% This function uses a modified version of emdos below, in order to choose
% the initial guess.


N = 4*1024;
[s s1 s2 s3] = gen_tests();


% General options
opt.maxmodes = 1;
opt.sporder = 8;
opt.alpha= 0.1;
opt.stop = 'f';
opt.orderder = [2 0 0]; % Selection derivee
opt.method = 'os';
opt.liss = [0];


p_alpha = [0.02:0.003:0.06 0.06:0.01:0.2];
%p_alpha = 0.04:0.01:0.05;

mesmet = {'integral mean','cubic spline','B-spline'};


err = zeros(length(mesmet),length(p_alpha));


for k=1:length(mesmet)
    for j=1:length(p_alpha)
        opt.alpha = p_alpha(j);
        [imf,ort,nder] = myemdos(mesmet{k},s,opt);
        %figure();plot(imf(1,:));
        err(k,j) = norm(imf(1,:)-s1);
    end
end

figure();
hold on;
plot(p_alpha,(1/N*err(1,:)),'ob-');
plot(p_alpha,(1/N*err(2,:)),'sr-');
plot(p_alpha,(1/N*err(3,:)),'*k-');

legend(mesmet{:});

xlabel('\alpha');
set(findall(0,'type','text'),'FontSize',15);
set(gca,'FontSize',15);
set(findall(0,'type','line'),'linewidth',2);
%set(gca,'Ylim',[0 0.025]);

end

function [imf,ort,nder] = myemdos(mamet,varargin)
%% EMDOS : function which computes the EMD by Optimisation on Splines.


[s,stop,alpha,maxmodes,t,method,orderder,liss,sporder] = init(varargin{:});


k = 1;
r=s;
imf = [];

%main loop : requires at least 3 extrema to proceed
while ~stop_emd(r) && (k < maxmodes+1 || maxmodes == 0)
    ar=r;
    switch(method)
        
        %% OS method
        case 'os'
            % Computation of m0
            [indexC tau m0 m0eval aux] = mybuild_m0(r,t,orderder(k),liss(k),sporder,mamet);
            
            nder(k) = aux;
            knots = aptknt(tau,sporder);
            L = length(indexC);
            m0 = m0';
            
            ti = t(indexC);
            D = matrice_der(knots,2,sporder);
            pass = matrice_pass(knots,ti,sporder);
            A = build_lambda(ti,sporder);
            A1 = (eye(L,L)+A)*pass;
        	b1 = alpha*abs((eye(L,L)-A)*(r(indexC)'-pass*m0));
            auxb1 = (eye(L,L)+A)*(r(indexC)');
            bb=vertcat(b1-auxb1,b1+auxb1);
            AA=vertcat(-A1,A1);
    
            % Computes x0 in the set of constraints
            x0 = pass\(s(indexC)');
            if sum(isnan(x0))>=eps
                x0 = [];
            end
    
            optis = optimset('MaxIter',10000,'LargeScale','off');
            [m,fval,exitflag] = quadprog(D,[],AA,bb,[],[],[],[],x0,optis);
    
        
            % Evaluation
            msp = spmak(knots,m');
            meval = fnval(msp,t);
            m = meval;
            r = r-m;
            
           
        case 'emd'
            noy = [1 2 1]/4;
            tmp=r;
            for j=1:liss(k)
                tmp = conv(tmp,noy,'valid');
            end
            
            % Prolongement
            [indmin,indmax] = extr(tmp);
            indmin=indmin+floor(liss(k));
            indmax=indmax+floor(liss(k));
            
            % Sifting
            stop_sift=0;
            aux=0;
            
            while ~stop_sift
                [tmin,tmax,mmin,mmax] = boundary_conditions(indmin,indmax,t,r,r,6);
                envmin = interp1(tmin,mmin,t,'spline');
                envmax = interp1(tmax,mmax,t,'spline');
                envmoy = (envmin+envmax)/2;
                nr = r-envmoy;
                
                
                switch(stop)
                    case 'f'
                        % Flandrin
                        amp = mean(abs(envmax-envmin))/2;
                        sx = abs(envmoy)./amp;
                        stop_sift = ~(mean(sx > alpha) > 0.05 | any(sx > 10*alpha));
                    case 'h'
                        % Huang
                        stop_sift = norm(nr-r)/(norm(r)+eps) < alpha;
                end
                
                if ~stop_sift
                    r=nr;
                    aux=aux+1;
                    [indmin,indmax] = extr(r);
                end
                nder(k)=aux;
            end
           
            
         case 'emdni'
             % Computation of m0
            [indmin indmax] = build_m02(r,t,orderder(k),liss(k));
            
            % Sifting
            stop_sift=0;
            aux=0;
            
            while ~stop_sift
                [tmin,tmax,mmin,mmax] = boundary_conditions(indmin,indmax,t,r,r,6);
                envmin = interp1(tmin,mmin,t,'spline');
                envmax = interp1(tmax,mmax,t,'spline');
                envmoy = (envmin+envmax)/2;
                nr = r-envmoy;
                
                
                switch(stop)
                    case 'f'
                        % Flandrin
                        amp = mean(abs(envmax-envmin))/2;
                        sx = abs(envmoy)./amp;
                        stop_sift = ~(mean(sx > alpha) > 0.05 | any(sx > 10*alpha));
                    case 'h'
                        % Huang
                        stop_sift = norm(nr-r)/(norm(r)+eps) < alpha;
                end
                
                if ~stop_sift
                    r=nr;
                    aux=aux+1;
                    %[indmin,indmax] = extr(r); % Sans lissage
                    [indmin indmax] = build_m02(r,t,orderder(k),liss(k)); % Avec lissage
                end
                nder(k)=aux;
            end
            

        otherwise
            error('Attention : wrong value of field ''method''');
    end

    imf(k,:) = r';
    r=ar-r;
    k=k+1;  
end
ort = io(s,imf);

% Residu
imf(k,:) = r';


end



function [indexC tau m0 m0eval aux] = mybuild_m0(r,t,ordersel,liss,sporder,mamet)
%% build_h0 : builds the approximation m0 of the local mean of the signal s
% INPUTS : 
%   s : signal
%   t : time
%   ordersel : manual selection of derivation order
%   spo : order of the splines
% OUTPUTS :
%   indexC : index of the extrema estimates
%   tau : points of the subdivision
%   m0 : estimation of the local mean : coefficients of the B-spline on tau
%   aux : effective order of serivation


%% Lissage du signal
noy = [1 2 1]/4;
tmp=r;
for j=1:liss
    tmp = conv(tmp,noy,'valid');
end
s=tmp;


%% Estimation of extrema of s
sder2 = madiff(s,t,2);
sder4 = madiff(s,t,4);
switch(ordersel)
    case -1
        seuil=2;
        % Automative selection
        [indmin, indmax, indzer] = extr(s,t);
        ls0 = length(indmin) + length(indmax);
        [indmin, indmax, indzer] = extr(sder2);
        ls2 = length(indmin) + length(indmax);
        [indmin, indmax, indzer] = extr(sder4);
        ls4 = length(indmin) + length(indmax);
        if ls2<=ls0+seuil
            [indmin, indmax, indzer] = extr(s,t);
            aux = 0;
        elseif ls4 <=ls2+seuil
            [indmin, indmax, indzer] = extr(sder2,t);
            aux = 2;
        else
            [indmin, indmax, indzer] = extr(sder4,t);
            aux = 4;
        end
        % Case discretization pb
        if max(ls0,max(ls2,ls4))>=length(s)/4
            warning('Pb de discretisation dans le calcul de la derivee');
            % Cas non C2
            if ls2<=ls0+seuil || ls2>=length(s)/4
                    [indmin, indmax, indzer] = extr(s,t);
                    aux = 0;
            else
                    [indmin, indmax, indzer] = extr(sder2,t);
                    aux = 2;
            end
        end

    case 0
        % Order 0
        [indmin, indmax, indzer] = extr(s);
        aux = 0;
    case 2
        % Order 2
        [indmin, indmax, indzer] = extr(sder2);
        aux = 2;
    case 4
        % Order 4
        [indmin, indmax, indzer] = extr(sder4);
        aux = 4;
    otherwise
        error('The ordersel field must be 0, 2 or 4');
end


%% Index of extremas
s=r;
indexC = sort([indmin indmax]);

% Si lissage
indexC = indexC+liss;
L = length(indexC);

xbar = zeros(1,L+1);
tbar = zeros(1,L+1);

switch(mamet)
    case 'integral mean'
        % Hong formula
        for i=1:L-1
            xbar(i+1) = mean(s(indexC(i)+1:indexC(i+1)-1));
            tbar(i+1) = floor(mean(t(indexC(i)+1:indexC(i+1)-1) .* ((s(indexC(i)+1:indexC(i+1)-1) -xbar(i+1)).^2))/...
            mean((s(indexC(i)+1:indexC(i+1)-1) -xbar(i+1)).^2));
        end

    case 'cubic spline'
        % Huang formula
        mamoy = 0.5*(interp1(t(indmin),r(indmin),t,'spline')+interp1(t(indmax),r(indmax),t,'spline'));
        indexbar = floor(0.5*(indexC(1:end-1)+indexC(2:end)));
        tbar(2:end-1) = t(indexbar);
        xbar(2:end-1) = mamoy(indexbar);
        
    case 'B-spline'
        % Huang B-spline formula
        valex = 0.25*(2*r(indexC(2:end-1))+r(indexC(1:end-2))+r(indexC(3:end)));
        masp = spapi(4,t(indexC(2:end-1)),valex);
        indexbar = floor(0.5*(indexC(1:end-1)+indexC(2:end)));
        mamoy = fnval(masp,t(indexbar));
        tbar(2:end-1) = t(indexbar);
        xbar(2:end-1) = mamoy;
        
end
            
            
            
        
        
        
% Symmetry
tbar(1) = 2*t(1)-tbar(2);
xbar(1) = xbar(2);
tbar(L+1) = 2*t(end)-tbar(L);
xbar(L+1) = xbar(L);




tau = tbar;
knots = aptknt(tau,sporder);
m0sp = spapi(knots,tau,xbar);
m0eval = fnval(m0sp,t);
m0 = m0sp.coefs;

end



