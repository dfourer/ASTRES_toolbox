function imellipsis = draw_ellipsis_rotated(xc,yc,a,b,theta,value,size)

N=size;

ellipsis=zeros(N,N);


%draw centered ellipsis
for x=floor(N/2)-a:floor(N/2)+a
   if (x >= 1 & x <= N)
       ymin = floor(N/2-b*sqrt(1-((x-N/2)/a)^2));
       ymax = floor(N/2+b*sqrt(1-((x-N/2)/a)^2));
       for y = ymin:ymax
           if (y >= 1 & y <= N)
               ellipsis(x,y)=value;
           end
       end
   end
end

%rotate ellipsis

ellipsis=imrotate(ellipsis,theta,'crop');

%translate ellipsis

imellipsis=zeros(N,N);
for i=1:N
    for j=1:N
        itranslate=floor(i+(xc-N/2));
        jtranslate=floor(j+(yc-N/2));
        if (itranslate >= 1 & itranslate <= N & jtranslate >= 1 & jtranslate <= N)
            imellipsis(itranslate,jtranslate)=ellipsis(i,j);
        end
    end
end
