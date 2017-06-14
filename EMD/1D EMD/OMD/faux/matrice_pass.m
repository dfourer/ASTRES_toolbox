function res = matrice_pass(sub,tau2,spo)
% matrice_pass : computes the matrix res such that, for any spline h(t) on a subdivision x and with
%   coefficients the vector H,
%       (M H)_j = h(x_i)
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr


% Nombre de splines
N = length(sub)-spo;
M = length(tau2);

res = zeros(M,N);

% Création matrice
for i=1:N
    % Création bspline j
    sp1 = spmak(sub(i:i+spo),1);
    res(:,i) = fnval(sp1,tau2');
    
end


end

