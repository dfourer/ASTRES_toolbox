%%method :  'Classic'
%           'OMD'
%           'VMD'
%           'Prox-EMD'
            

function modes = emd_1D_main(x,nc,method,varargin)

switch nargin
    case 1
        nc = 0;
        method = 'Classic';
    case 2
        method = 'Classic';
end
    
if strcmp(method,'Classic')
    if nc == 0
        modes = emd(x);
    else modes = emd(x,'maxmodes',nc-1);
    end
end

if strcmp(method,'OMD')
    if nc == 0
        modes = emdos(x);
    else modes = emdos(x,'mexmodes',nc-1);
    end 
end

if strcmp(method,'VMD')
     if nc == 0
        x_emd  = emd(x');
        nc = size(x_emd,1);
     end
     % Default values
     alpha = 1;
     tau = 1;
     DC = 0;
     init = 0;
     switch length(varargin)
         case 1
             alpha = varargin{1};
         case 2
             alpha = varargin{1};
             tau = varargin{2};
         case 3
             alpha = varargin{1};
             tau = varargin{2};
             DC = varargin{3};
         case 4
             alpha = varargin{1};
             tau = varargin{2};
             DC = varargin{3};
             init = varargin{4};
     end
     tol = 1e-6;
     modes = VMD(x, alpha, tau, nc, DC, init, tol);
end

if strcmp(method,'Prox-EMD')
    if nc == 0
        x_emd  = emd(x');
        nc = size(x_emd,1);
    end
    data.K = nc;
    algo=[];
     switch length(varargin)
         case 1
             algo.epsa = varargin{1};
         case 2
             algo.epsa = varargin{1};
             algo.epsd = varargin{2};
     end
    modes = prox_emd(x, data, algo);
end

end