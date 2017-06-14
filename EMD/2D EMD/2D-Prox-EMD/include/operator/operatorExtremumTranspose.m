function yt = operatorExtremumTranspose(y,A)
[nx,ny]=size(y);
n=nx*ny;

y=reshape(y,1,n);

yt = A'*y';
%yt = yt';
yt=reshape(yt,nx,ny);