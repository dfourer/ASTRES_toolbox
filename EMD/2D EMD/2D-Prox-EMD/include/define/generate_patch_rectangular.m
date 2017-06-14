function patch = generate_patch_rectangular(x,y,lx,ly,intensity,size)

N=size;

patch = zeros(N,N);

xmin=max(x-lx/2,1);
xmax=min(x+lx/2,N);

ymin=max(y-ly/2,1);
ymax=min(y+ly/2,N);

patch(xmin:xmax,ymin:ymax)=intensity;

end