function results = run_prox_emd_2D(data,type,lambda)
%
% The function run_prox_emd_2D provides a 2D mode decomposition of the
% signal data.x of size (data.n,data.m). The decomposition is performed
% according to the work "2D Prony-Huang Transform: A New Tool for 2-D
% Spectral Analysis, J. Schmitt, N. Pustelnik, P. Borgnat, P. Fladrin, L.
% Condat, IEEE Transactions on Image Processing, 2014".
%
% The algorithm decomposes the data x into a trend component a (results.tend)
% and an IMF component d (results.imf).
% It consists to solve
% For k=1 : K-1
%     (a_k,d_k) = \arg\min_{a,d} ||x-a-d||_2^2 + lambda_geometry.||A d||_{2,1} + lambda.texture.||D d ||_1
%     x = a_k
% End
% The parameters involved in the criterion can be tuned in order to provide
% good estimates.
%
%%%%%%%%%%%%INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% data                  : data to decompose
% data.x                : image
% data.n, data.m        : size of image
%
% type                  : choice of the type of EMD
%                         4 approaches are available :
%                         'G2D'          Genuine 2D
%                         'G2Ddir'       Directional genuine 2D
%                         'P2DLC'        Pseudo 2D on lines and columns
%                         'P2DLCD'       Pseudo 2D on lines, columns, diagonals and anti-diagonals
%                         For differences between these approaches, see 2D Prony-Huang Transform: A New Tool for 2-D
%                         Spectral Analysis, J. Schmitt, N. Pustelnik, P. Borgnat, P. Fladrin, L.
%                         Condat, IEEE Transactions on Image Processing, 2014"
%
% lambda                 : parameters of algorithm
% lambda.K               : number of IMFs
% lambda.lambda_geometry : cell of size K containing the regularization
%                          parameters for each trend
% lambda.lambda_texture(_l, _c, _d1, d2) : cell of size K containing the 
%                          regularization parameters for each IMF
%
%%%%%%%%%%%%OUTPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% results                : results of algorithm
% results.trend          : trend of order K
% results.imf            : cell of size K containing the K IMFs
% results.crit           : cell of size K containing the variation of the
%                          convergence criterion 
%                          ||x-a-d||_2^2 + lambda_geometry.||A d||_{2,1} + lambda.texture.||D d ||_1
%
%%%%%%%%%%%%SUBFUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%
%
% define_method_emd_2D  : for each mode, define parameters, functions, operators used in the iterative algorithm 
% 
% The 3 other subfunctions perform the decomposition under the 4 different approaches :
% prox_emd_2D_G2D       : Genuine 2D and directional genuine 2D
% prox_emd_2D_P2DLC     : Pseudo 2D on lines and columns
% prox_emd_2D_P2DLCD    : Pseudo 2D on lines, columns, diagonals and anti-diagonals

for k=1:lambda.K
    [method,param,dec,func,grad_prox] = define_method_emd_2D(data,type,lambda,k);

    if strcmp(method.filter2, 'P2DLC')
        resimf = prox_emd_2D_P2DLC(data, dec, func, grad_prox, param, k);
        elseif strcmp(method.filter2, 'P2DLCD')
        resimf = prox_emd_2D_P2DLCD(data, dec, func, grad_prox, param, k);
        else resimf = prox_emd_2D_G2D(data, dec, func, grad_prox, param, k);   
    end
    results.imf{k} = resimf.imf;
    results.crit{k} = resimf.crit;
    data.x = resimf.trend;
end
results.trend = resimf.trend;

end