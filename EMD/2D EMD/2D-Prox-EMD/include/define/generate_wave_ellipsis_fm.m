function wave = generate_wave_ellipsis_fm (amp,fp,angle,fm,m,xc,yc,a,b,theta,size)

N=size;
theta = theta*180/pi;
ellipsis = draw_ellipsis_rotated(xc,yc,a,b,theta,1,size);

wave = zeros(N,N);
for i=1:N
    for j=1:N
        wave(i,j) = amp*cos(2*pi*fp*(cos(angle)*i+sin(angle)*j)/(N) + m*sin((2*pi*fm*(cos(angle)*i+sin(angle)*j))/N))*ellipsis(i,j);
    end
end
%figure;
%imagesc(wave);colormap(gray);

end