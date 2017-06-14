function yt = operatorExtremum(y,A)
[nx,ny]=size(y);
n=nx*ny;

y=reshape(y,n,1);
yt=A*y;
%yt = A*y';
%yt = yt';
yt=reshape(yt,nx,ny);