%extracts the indices of extrema
% Writtezn by Gabriel Rilling
function [indmin, indmax, indzer] = extr(x,t)

if(nargin==1)
  t=1:length(x);
end

m = length(x);

if nargout > 2
  x1=x(1:m-1);
  x2=x(2:m);
  indzer = find(x1.*x2<0);

  if any(x == 0)
    iz = find( x==0 );
    indz = [];
    if any(diff(iz)==1)
      zer = x == 0;
      dz = diff([0 zer 0]);
      debz = find(dz == 1);
      finz = find(dz == -1)-1;
      indz = round((debz+finz)/2);
    else
      indz = iz;
    end
    indzer = sort([indzer indz]);
  end
end

d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;


% when two or more successive points have the same value we consider only one extremum in the middle of the constant area
% (only works if the signal is uniformly sampled)

if any(d==0)

  imax = [];
  imin = [];

  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if length(debs) > 0
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)];
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)];
        end
      end
    end
  end

  if length(imax) > 0
    indmax = sort([indmax imax]);
  end

  if length(imin) > 0
    indmin = sort([indmin imin]);
  end

end
end
