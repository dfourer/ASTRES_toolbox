function imellipsis = draw_ellipsis(xc,yc,a,b,value,size)

N=size;

imellipsis=zeros(N,N);

for i=xc-a:xc+a
    if (i > 1 & i < N)
        jmin = floor(yc-b*sqrt(1-((i-xc)/a)^2));
        jmax = floor(yc+b*sqrt(1-((i-xc)/a)^2));
        for j = jmin:jmax
            if (j > 1 & j < N)
                imellipsis(i,j)=value;
            end
        end
    end
end

