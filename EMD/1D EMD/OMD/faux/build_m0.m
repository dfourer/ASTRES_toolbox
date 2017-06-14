function [indexC tau m0 m0eval aux] = build_m0(r,t,ordersel,sporder)
%% build_h0 : builds the approximation m0 of the local mean of the signal s
% INPUTS : 
%   r : signal
%   t : time
%   ordersel : manual selection of derivation order
%   sporder : order of the splines
% OUTPUTS :
%   indexC : index of the extrema estimates
%   tau : points of the subdivision
%   m0 : estimation of the local mean : coefficients of the B-spline on tau
%   m0eval : the B-spline evaluated on t
%   aux : effective order of derivation
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr

s=r;

%% Estimation of extrema of s
sder2 = madiff(s,t,2);
sder4 = madiff(s,t,4);
switch(ordersel)
    case -1
        seuil=8;
        % Automative selection (see [1])
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
            % Case not C2
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

L = length(indexC);

% Hong formula
xbar = zeros(1,L+1);
tbar = zeros(1,L+1);
for i=1:L-1
	xbar(i+1) = mean(s(indexC(i):indexC(i+1)));
    tbar(i+1) = mean(t(indexC(i):indexC(i+1)) .* ((s(indexC(i):indexC(i+1)) -xbar(i+1)).^2))/...
    mean((s(indexC(i):indexC(i+1)) -xbar(i+1)).^2);
end 

% Symmetry
tbar(1) = 2*t(1)-tbar(2);
xbar(1) = xbar(2);
tbar(L+1) = 2*t(end)-tbar(L);
xbar(L+1) = xbar(L);


tau = tbar;
knots = aptknt(tau,sporder);
m0sp = spapi(knots,tau,[xbar(1) 0 xbar(3:end-1) 0]);
m0eval = fnval(m0sp,t);
m0 = m0sp.coefs;

end

