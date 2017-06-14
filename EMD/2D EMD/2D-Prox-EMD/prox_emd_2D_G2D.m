function resimf = prox_emd_2D_G2D(data, dec, func, grad_prox, param, k)

x1      = zeros(size(data.x));
x2      = zeros(size(data.x));
y11      = dec.T1(zeros(size(data.x)));
y12      = dec.T2(zeros(size(data.x)));

sigma   = param.sigma;
tau     = param.tau;
rho     = param.rho;

lambda_geometry = param.lambda_geometry;
lambda_texture = param.lambda_texture;

crit  = zeros(1,param.Lmax);

param.naff=1000;

for n = 1:param.Lmax
    
    x1tilde = x1 - tau*(grad_prox.g(x1,x2) + dec.T1_adj(y11));
    x2tilde = x2 - tau*(grad_prox.g(x1,x2) + dec.T2_adj(y12));
    
    y11tilde = y11 + sigma*(dec.T1(2*x1tilde - x1));
    y11tilde = y11tilde - sigma*grad_prox.f1(y11tilde/sigma,lambda_geometry/sigma);
    
    y12tilde = y12 + sigma*(dec.T2(2*x2tilde - x2));
    y12tilde = y12tilde - sigma*grad_prox.f2(y12tilde/sigma,lambda_texture/sigma);
 
    
    x1  = rho*x1tilde  + (1-rho)*x1;
    x2  = rho*x2tilde  + (1-rho)*x2;
    y11 = rho*y11tilde + (1-rho)*y11; 
    y12 = rho*y12tilde + (1-rho)*y12; 
    
    
    %evaluate the criterion
    crit(n)  = func.g(x1,x2) + func.f1(x1,lambda_geometry) + func.f2(x2,lambda_texture);
    
    %stop criterion
    if n>2
       if abs(crit(n)- crit(n-1))/abs(crit(n) +eps)<1e-6
           break;
       end
    end
    
    
end
resimf.trend      = x1;
resimf.imf        = x2;
resimf.residual   = data.x - x1 - x2;
resimf.crit  = crit(1:n)



end