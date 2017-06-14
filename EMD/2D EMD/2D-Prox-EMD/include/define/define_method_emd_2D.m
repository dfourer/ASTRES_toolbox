function [method,param,dec,func,grad_prox] = define_method_emd_2D(data,type,lambda,k)

method.filter1        = 'TV';
method.filter2        = type;

if any([strcmp(method.filter2,'G2D');strcmp(method.filter2,'G2Ddir')])
    param.Lmax            = 10000;
else param.Lmax = 3000;
end

if strcmp(method.filter1,'TV')  
    param.lambda_geometry = lambda.lambda_geometry{k};
    h = [1/2 -1/2; 0 0];
    H = psf2otf(h,[data.n,data.m]);
    v = h';
    V = psf2otf(v,[data.n,data.m]);  
    T = zeros(data.n,data.m,2);
    T(:,:,1)=H;
    T(:,:,2)=V;
    dec.T1       = @(y) operatorTV(y,T);
    dec.T1_adj   = @(y) operatorTVTranspose(y,T);
end

if strcmp(method.filter2,'G2D')
    param.lambda_texture = lambda.lambda_texture{k};
    [D,normD,~,~,~] = defineParameterMax_2D(data.x); 
    dec.T2       = @(y) operatorExtremum(y,D);
    dec.T2_adj   = @(y) operatorExtremumTranspose(y,D);
    param.normA = 2.4165;
    elseif strcmp(method.filter2, 'G2Ddir')
    param.lambda_texture = lambda.lambda_texture{k};
    [D,normD,~,~,~] = defineParameterMax_2D_directional(data.x); 
    dec.T2       = @(y) operatorExtremum(y,D);
    dec.T2_adj   = @(y) operatorExtremumTranspose(y,D);
    param.normA = 2.4165;

    
elseif strcmp(method.filter2, 'P2DLC')
    param.lambda_texture_l = lambda.lambda_texture_l{k};
    param.lambda_texture_c = lambda.lambda_texture_c{k};
    [Dl,normDl] = defineParameterMax_line(data.x);
    [Dc,normDc] = defineParameterMax_column(data.x);
    dec.T2      = @(y) operatorExtremum(y,Dl);
    dec.T2_adj   = @(y) operatorExtremumTranspose(y,Dl);
    dec.T3       = @(y) operatorExtremum(y,Dc);
    dec.T3_adj   = @(y) operatorExtremumTranspose(y,Dc);
    param.normA = 2.8618;
elseif strcmp(method.filter2, 'P2DLCD')
    param.lambda_texture_l = lambda.lambda_texture_l{k};
    param.lambda_texture_c = lambda.lambda_texture_c{k};
    param.lambda_texture_d1 = lambda.lambda_texture_d1{k};
    param.lambda_texture_d2 = lambda.lambda_texture_d2{k};
    [Dl,normDl] = defineParameterMax_line(data.x);
    [Dc,normDc] = defineParameterMax_column(data.x);
    [Dd1,normDd1] = defineParameterMax_diagonal1(data.x);
    [Dd2,normDd2] = defineParameterMax_diagonal2(data.x);
    dec.T2      = @(y) operatorExtremum(y,Dl);
    dec.T2_adj   = @(y) operatorExtremumTranspose(y,Dl);
    dec.T3       = @(y) operatorExtremum(y,Dc);
    dec.T3_adj   = @(y) operatorExtremumTranspose(y,Dc);
    dec.T4       = @(y) operatorExtremum(y,Dd1);
    dec.T4_adj   = @(y) operatorExtremumTranspose(y,Dd1);
    dec.T5       = @(y) operatorExtremum(y,Dd2);
    dec.T5_adj   = @(y) operatorExtremumTranspose(y,Dd2);
    param.normA = 4.1416;
end

func.g       = @(y1,y2) norm(data.x - y1-y2, 'fro')^2;  
grad_prox.g  = @(y1,y2) 2*(y1 + y2 - data.x); 

func.f1       = @(y, threshold) threshold * sum(reshape(sqrt(sum(abs(dec.T1(y)).^2,3)), 1, data.n*data.m));  
grad_prox.f1  = @(y, threshold) prox_L12(y, threshold);  
func.f2      = @(y, threshold) threshold * sum(abs(reshape(dec.T2(y), 1, numel(dec.T2(data.x)))));  
grad_prox.f2 = @(y, threshold) prox_L1(y, threshold);



if strcmp(method.filter2, 'P2DLC')
    func.f3      = @(y, threshold) threshold * sum(abs(reshape(dec.T3(y), 1, numel(dec.T3(data.x)))));    
    grad_prox.f3 = @(y, threshold) prox_L1(y, threshold);
end

if strcmp(method.filter2, 'P2DLCD')
    func.f3      = @(y, threshold) threshold * sum(abs(reshape(dec.T3(y), 1, numel(dec.T3(data.x)))));    
    grad_prox.f3 = @(y, threshold) prox_L1(y, threshold);
    func.f4      = @(y, threshold) threshold * sum(abs(reshape(dec.T4(y), 1, numel(dec.T4(data.x)))));    
    grad_prox.f4 = @(y, threshold) prox_L1(y, threshold);
    func.f5      = @(y, threshold) threshold * sum(abs(reshape(dec.T5(y), 1, numel(dec.T5(data.x)))));    
    grad_prox.f5 = @(y, threshold) prox_L1(y, threshold);
end
    
param.sigma = 1;
param.tau = 0.99/(param.sigma*(param.normA) + 2 );
param.rho = 1;

end
