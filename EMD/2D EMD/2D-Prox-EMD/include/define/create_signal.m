function data=create_signal(type)

if strcmp(type,'example1')
    wave1 = generate_wave(10,2,-pi/6,256);
    wave2 = generate_wave(10,10,pi/4,256);
    data.x = wave1 + wave2;
    [nx,ny] = size(data.x);
    data.n = nx;
    data.m = ny;
    clear nx ny;
end

if strcmp(type,'example2')
    wave1 = generate_wave(10,0.5,pi/2,256);
    wave2 = generate_wave(10,5,pi/4,256);
    %wave3 = generate_wave(10,40,-pi/4,256);
    wave3 = generate_wave(10,40,pi/3,256);
    data.x = wave1 + wave2 + wave3;
    [nx,ny] = size(data.x);
    data.n = nx;
    data.m = ny;
    clear nx ny;
end

if strcmp(type,'example3')
    wave1 = generate_wave_ellipsis_fm(20,60,pi/3,10,1,200,250,170,75,pi/3,512);
    wave2 = generate_wave_ellipsis_fm(20,120,pi/6,20,1,300,200,170,75,-pi/6,512);
    patch1 = generate_patch_rectangular(400,250,75,200,10,512);
    patch2 = draw_ellipsis(100,100,40,30,10,512);
    data.x = wave1+wave2+patch1+patch2;
    [nx,ny] = size(data.x);
    data.n = nx;
    data.m = ny;
    clear nx ny;
end