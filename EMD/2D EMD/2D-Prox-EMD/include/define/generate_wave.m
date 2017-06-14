function wave = generate_wave (amp,freq,angle,size)

N=size;

wave = zeros(N,N);
for i=1:N
    for j=1:N
        wave(i,j) = amp*cos(2*pi*freq*(cos(angle)*i+sin(angle)*j)/(N));
    end
end
%imagesc(wave);colormap(gray);

end