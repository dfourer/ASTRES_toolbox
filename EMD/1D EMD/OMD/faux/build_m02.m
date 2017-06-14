function [indmin indmax] = build_m02(r,t,ordersel,liss)
%% build_h0 : builds the approximation m0 of the local mean of the signal s
% INPUTS : 
%   s : signal
%   t : time
%   ordersel : manual selection of derivation order
%   spo : order of the splines
% OUTPUTS :
%   indexC : index of the extrema estimates
%   tau : points of the subdivision
%   m0 : estimation of the local mean : coefficients of the B-spline on tau
%   aux : effective order of serivation
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr


%% Lissage du signal
noy = [1 2 1]/4;
tmp=r;
for j=1:liss
    tmp = conv(tmp,noy,'valid');
end
s=tmp;


%% Estimation of extrema of s
sder2 = madiff(s,t,2);
sder4 = madiff(s,t,4);
switch(ordersel)
    case -1
        seuil=8;
        % Automative selection
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
            % Cas non C2
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

indmin = indmin+liss;
indmax = indmax+liss;

% Si bruit
if 0%length(indmin)+length(indmax)>length(s)/2
    warning('Using sure SHRINKAGE');
    sd = wden(s,'heursure','s','one',3,'sym8');
    [indmin indmax] = extr(s-sd);
    %indexC = sort([indmin indmax]);
end

end

