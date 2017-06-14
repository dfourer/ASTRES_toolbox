function results = prox_emd(x, data, algo)
%
%
% The function "prox_emd.m" provides a mode decomposition of the signal x
% of length N. The decomposition is performed according to the work 
% "N. Pustelnik, P. Borgnat, P. Flandrin, Empirical Mode Decomposition 
% revisited by multicomponent non smooth convex optimization, Signal 
% Processing, Vol. 102, pp. 313--331, Sept. 2014."
%
% It consists to solve
% For k=1 : K-1
%     (a_k,d_k) = \arg\min_{a,d} ||x-a-d||^p + c1.||A a||^q + c2.||D d ||^r
%     x = a_k
% End
% The parameters involved in the criterion can be tuned in order to provide
% good estimates.
%
%
% data.K           : number of components K
%                    Default:  value computed with usual EMD
%
% algo.nbiter      : number of iteration for the resolution of the 
%                    minimization problem for each k 
%                    Example: 
%                       >> algo.nbiter = [3000, 10000];
%                    Default:
%                       >> algo.nbiter = 3000*ones(1,data.K);
%
% algo.epsa        : Value of c1
%                    Example: 
%                       >> algo.epsa = [19,0.1];
%                    Default: value of ||A a||^q when a denotes a solution 
%                             from usual EMD
%
% algo.epsd        : Value of c2
%                    Example: 
%                       >> algo.epsd = [0.05,0.01];
%                    Default: value of ||D d||^r when a denotes a solution 
%                             from usual EMD
%
% algo.dataFid      : Value of p ('l1' or 'l2');
%                    Example: 
%                       >> algo.dataFid{1}  = 'l2';
%                       >> algo.dataFid{2}  = 'l1';
%                    Default: 
%                       >> algo.dataFid  = 'l2';
%
% algo.regTre      : Value of q ('l1' or 'l2');
%                    Example: 
%                       >> algo.regTre{1}  = 'l2';
%                       >> algo.regTre{2}  = 'l1';
%                    Default: 
%                       >> algo.regTre  = 'l2';
%
% algo.regFlu      : Value of r ('l1' or 'l2');
%                    Example: 
%                       >> algo.regFlu{1}  = 'l2';
%                       >> algo.regFlu{2}  = 'l1';
%                    Default: 
%                       >> algo.regFlu  = 'l1';
%
% algo.filt        : Value of A ('laplacian' or 'gradient');
%                    Example: 
%                       >> algo.regFlu{1}  = 'laplacian';
%                       >> algo.regFlu{2}  = 'gradient';
%                    Default: 
%                       >> algo.regFlu  = 'laplacian';
%
% algo.ext         : Value of D ('d1' or 'd2');
%                    Example: 
%                       >> algo.ext{1}  = 'd1';
%                       >> algo.ext{2}  = 'd1';
%                    Default: 
%                       >> algo.ext  = 'd1';
%
%
% results          : returns the components after decomposition . It consists in a
%                    matrix of size KxN where K denotes the number of extracted
%                    components. [d_1;d_2;...d_{K-1};a_{K-1}]
% 
%
% Basic use of prox_emd :
% >> results     = prox_emd(x);

if nargin<3
    algo = [];
else
    algo
end
data.x = x;
data.n = length(data.x);
x_emd  = emd(data.x');


if ~isfield(data,'K')
    data.K = size(x_emd,1);
else
    x_emd = emd(data.x','maxmodes',data.K-1);
end
results = zeros(data.K,size(data.x,2));

if ~isfield(algo,'dataFid')
    for k=1:data.K-1
        algo.dataFid{k} = 'l2';
    end
end
if ~isfield(algo,'regTre')
     for k=1:data.K-1
        algo.regTre{k} = 'l2';
     end
end
if ~isfield(algo,'filt')
     for k=1:data.K-1
        algo.filt{k}= 'laplacian';
     end
end
if ~isfield(algo,'regFlu')
     for k=1:data.K-1
        algo.regFlu{k} = 'l1';
     end
end
if ~isfield(algo,'ext')
     for k=1:data.K-1
        algo.ext{k} = 'd1';
     end
end
if ~isfield(algo,'nbiter')
    algo.nbiter = 3000*ones(1,data.K-1);
end


if ~isfield(algo,'epsa')
    algo.epsa = NaN*ones(1,data.K-1);
end
if ~isfield(algo,'epsd')
    algo.epsd = NaN*ones(1,data.K-1);
end


for k=1:data.K-1    
    algo.dataFidelity   = algo.dataFid{k};
    algo.regTrend       = algo.regTre{k};
    algo.filter         = algo.filt{k};    
    algo.regFluct       = algo.regFlu{k};    
    algo.extr           = algo.ext{k};
    
    emd_define_linear_operators;
    emd_define_parameters_algo;

    
    if isnan(algo.epsa(k))
        if strcmp(algo.regTrend,'l2')
            algo.epsa(k)  = norm(f1.L_a(x_emd(k,:)))^2/100; 
        elseif strcmp(algo.regTrend,'l1')
            algo.epsa(k)  = sum(abs(f1.L_a(x_emd(k,:))))/100;
        end
    end
    if isnan(algo.epsd(k))
        if strcmp(algo.regFluct,'l2')
            algo.epsd(k)  = norm(f1.L_d(sum(x_emd(k+1:end,:),1),D))^2/10; 
        elseif strcmp(algo.regFluct,'l1')
            algo.epsd(k)  = sum(abs(f1.L_d(sum(x_emd(k+1:end,:),1),D)))/10;
        end
    end

    algo.eps_a          = algo.epsa(k);
    algo.eps_d          = algo.epsd(k);
    algo.niter          = algo.nbiter(k);
    algo.theta          = max( sqrt(2),max(1,normA));
    algo.epsilon        = 0.90/(algo.theta+1);
    algo.gamma          = 1*(1-algo.epsilon)/algo.theta;
    emd_define_criteria;
    emd_define_proximity_operators;
    t                   = emd_algorithm(data,algo,f1,f2,D);
    data.x              = t(2,:);
    results(k,:)        = t(1,:);
    fprintf('>> IMF-Trend of order %d computed\n',k);
end
results(data.K,:)       = t(2,:);


