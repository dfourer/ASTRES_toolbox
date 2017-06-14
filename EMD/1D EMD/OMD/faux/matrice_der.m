function res = matrice_der(sub,order,spo)
%% MATRICE_DER : internal function of EMDOS. Computes the "differentiation matrix" for b-splines:
%       res(i,j)=int B^(d)_i(t) B^(d)_j(t) dt
%       where B^(d)_j is the d-derivative of the j-th B-spline
% INPUTS:
%   sub : knots subdivision
%   order : differentiation order
%   spo : spline order
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr

% Nombre de splines
N = length(sub)-spo;

res = zeros(N,N);

for i=1:N
    % Création bspline i
	sp1 = spmak(sub(i:i+spo),1);
    sp1 = fnder(sp1,order);
    
    for j=i:min(i+spo-1,N)
        % Création spline 2
        sp2 = spmak(sub(j:j+spo),1);
        sp2 = fnder(sp2,order);
        aux = fncmb(sp1,'*',sp2);
        res(i,j) = diff(fnval(fnint(aux),[sub(i) sub(i+spo)]))/N^order;
    end
end

% Complétion

% Triangle inférieur (par symétrie)
for i=2:N
    for j=max(1,i-spo-1):i-1
        res(i,j)=res(j,i);
    end
end

end
