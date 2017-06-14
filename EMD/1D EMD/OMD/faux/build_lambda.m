function C = build_lambda(ti,k)
%% build_lambda : internal function of EMDOS. Computes the matrix C st:
%   C*h(xi) = lambda_i
%  where h(xi) is the vector of the spline evaluated at the extrema xi, and
%  lambda_i is the symmetric points defined in [1].%
%
% Thomas Oberlin
% 12.2011
% thomas.oberlin@imag.fr

n=k-1;
r = k/2;


ti1 = ti(1:2:end);
ti2 = ti(2:2:end);


C1 = spcol([ti1(1)*ones(1,r) ti1 ti1(end)*ones(1,r)],k,ti2);
C2 = spcol([ti2(1)*ones(1,r) ti2 ti2(end)*ones(1,r)],k,ti1);
% "Mélange" des matrices C1 et C2.
C = zeros(size(C1)+size(C2));
for i=r:size(C1,1)-r
	% On commence par C1
    C(2*i,2*(i-r)+1:2:2*(i+r)) = C1(i,i-r+1:i+r);
    C(2*i+1,2*i+1-n:2:2*i+1+n) = C2(i,i-r+1:i+r);
end
% On complète (ordre 1)
for i=1:2*r-1
	C(i,i+1)=1;
end
for i=size(C,1)-2*r+2:size(C,1)
	C(i,i-1)=1;
end
    


end