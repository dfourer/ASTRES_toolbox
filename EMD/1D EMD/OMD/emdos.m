function [imf,ort,nder] = emdos(varargin)
%% EMDOS : Computes the EMD by Optimisation on Splines, implements the method described in [1]. Uses either the standard EMD [2,3] or the OS algorithm [1].
% 
% SYNTAX
%   imf = emdos(s)
%   imf = emdos(s,...,'Option_name',Option_value,...)
%   imf = emdos(s,opt)
%   [imf,ort,order_der] = emdos(...)
%
% INPUTS
%   s : one-dimensional real signal
%   opt : field of optional parameters. They are listed below together with
%   the default value
%       method      <'os'>  : 'emd' or 'emdni' or 'os'
%       stop        <'f'>   : Kind of stopping criterion : 'h' for the Huang one [2] and
%                               'f' for Flandrin and Rilling's one [3].
%       alpha       <0.05>  : Stopping criterion parameter [3]. 
%       maxmodes    <10>    : Maximum number of IMFs
%       t                   : Time vector t
%       orderder    <-1>    : Differentiation order for each mode : either a vector or a scalar.
%                               -1 : automatic procedure described in [1]
%                               0, 2 or 4 : uses the 0th, 2th or 4th
%                                   derivative respectively.
%       sporder     <8>     : Order of the splines used in the os method.
%                               sporder>=6 and <=12, must be even.
%
% OUTPUTS
%   imf     : (number of IMFs +1) X length(s) array which contains each IMFs
%               and a residue
%   ort     : orthogonality index defined in [3,1]
%   nder    : derivation order used for each IMF (array of length the
%               number of IMFs).
% EXAMPLE
%   n = 1024;
%   t = linspace(0,1,n);
%   s = cos(60*pi*(t+1).^2) + (1+exp(-18*(t-0.5).^2)).*sin(50*pi*t);
%   opt.maxmodes = 1;
%   opt.sporder = 8;
%   opt.alpha= 0.05;
%   opt.stop = 'f';
%   opt.orderder = [2 0]
%   opt.method = 'os';
%   opt.t = t;
%   [imf ort nder] = emdos(s,opt);
%   plot_imf(imf,t,'');
% 
% REFERENCES
%   [1] T. Oberlin, S. Meignen and V. Perrier, “An Alternative Formulation
%       for the Empirical Mode Decomposition”, IEEE Trans. on Signal
%       Prcessing, to appear.
%   [2] N. Huang, Z. Shen, S. Long, M. Wu, H. Shih, Q. Zheng, N. Yen, C. Tung, and H. Liu, “The empirical mode decomposition
%       and the Hilbert spectrum for nonlinear and non-stationary time series analysis,” Proceedings of the Royal Society :
%       Mathematical, Physical and Engineering Sciences, vol. 454, no. 1971, pp. 903–995, 1998.
%   [3] G. Rilling, P. Flandrin, and P. Gonc¸alv`es, “On empirical mode decomposition and its algorithms,” in IEEE-EURASIP
%       workshop on nonlinear signal and image processing NSIP-03, Grado (I), 2003.
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr





% Gets the parameter
[s,stop,alpha,maxmodes,t,method,orderder,liss,sporder] = init(varargin{:});

k = 1;
r=s;
imf = [];

%main loop : requires at least 3 extrema to proceed
while ~stop_emd(r) && (k < maxmodes+1 || maxmodes == 0)
    ar=r;
    switch(method)
        
        % OS method
        case 'os'
            % Computation of m0
            [indexC tau m0 m0eval aux] = build_m0(r,t,orderder(k),sporder);
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
    
            optis = optimset('MaxIter',20000,'LargeScale','off');
            [m,fval,exitflag] = quadprog(D,[],AA,bb,[],[],[],[],x0,optis);
    
            % Evaluation
            msp = spmak(knots,m');
            meval = fnval(msp,t);
            m = meval;
            r = r-m;
            
               
        case 'emd'
            tmp = r;
            % padding
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