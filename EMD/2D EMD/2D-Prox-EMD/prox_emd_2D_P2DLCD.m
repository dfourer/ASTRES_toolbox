function resimf = prox_emd_2D_P2DLCD(data, dec, func, grad_prox, param, k)

x1      = zeros(size(data.x));
x2      = zeros(size(data.x));
y11      = dec.T1(zeros(size(data.x)));
y12      = dec.T2(zeros(size(data.x)));
y13      = dec.T3(zeros(size(data.x)));
y14      = dec.T4(zeros(size(data.x)));
y15      = dec.T5(zeros(size(data.x)));


sigma = param.sigma;
tau = param.tau;
rho = param.rho;

lambda_geometry = param.lambda_geometry;
lambda_texture_l = param.lambda_texture_l;
lambda_texture_c = param.lambda_texture_c;
lambda_texture_d1 = param.lambda_texture_d1;
lambda_texture_d2 = param.lambda_texture_d2;


crit = zeros(1,param.Lmax);

for n = 1:param.Lmax
       
    x1tilde = x1 - tau*(grad_prox.g(x1,x2) + dec.T1_adj(y11)) ;
    x2tilde = x2 - tau*(grad_prox.g(x1,x2) + dec.T2_adj(y12) + dec.T3_adj(y13) + dec.T4_adj(y14) + dec.T5_adj(y15)); 
    
    y11tilde = y11 + sigma*(dec.T1(2*x1tilde - x1));
    y11tilde = y11tilde - sigma*grad_prox.f1(y11tilde/sigma,lambda_geometry/sigma);

    y12tilde = y12 + sigma*(dec.T2(2*x2tilde - x2));
    y12tilde = y12tilde - sigma*grad_prox.f2(y12tilde/sigma,lambda_texture_l/sigma);
    
    y13tilde = y13 + sigma*(dec.T3(2*x2tilde - x2));
    y13tilde = y13tilde - sigma*grad_prox.f3(y13tilde/sigma,lambda_texture_c/sigma);
    
    y14tilde = y14 + sigma*(dec.T4(2*x2tilde - x2));
    y14tilde = y14tilde - sigma*grad_prox.f4(y14tilde/sigma,lambda_texture_d1/sigma);
    
    y15tilde = y15 + sigma*(dec.T5(2*x2tilde - x2));
    y15tilde = y15tilde - sigma*grad_prox.f5(y15tilde/sigma,lambda_texture_d2/sigma);
    
    x1  = rho*x1tilde  + (1-rho)*x1;
    x2  = rho*x2tilde  + (1-rho)*x2;
    y11 = rho*y11tilde + (1-rho)*y11; 
    y12 = rho*y12tilde + (1-rho)*y12; 
    y13 = rho*y13tilde + (1-rho)*y13;
    y14 = rho*y14tilde + (1-rho)*y14; 
    y15 = rho*y15tilde + (1-rho)*y15; 


    
    %evaluate the criterion
    crit(n)   = func.g(x1,x2) + func.f1(x1,lambda_geometry) + func.f2(x2,lambda_texture_l) + func.f3(x2,lambda_texture_c) + func.f4(x2,lambda_texture_d1) + func.f5(x2,lambda_texture_d2);

    %stop criterion
    if n>2
        if abs(crit(n)- crit(n-1))/abs(crit(n) +eps)<1e-6
            break;
        end
     end

end

resimf.trend    = x1;
resimf.imf      = x2;
resimf.residual = data.x - x1 - x2;
resimf.crit  = crit(1:n);

end