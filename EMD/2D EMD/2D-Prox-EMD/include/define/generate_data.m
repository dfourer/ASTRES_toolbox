wave1=generate_wave_ellipsis_fm(20,60,pi/3,10,1,200,250,170,75,pi/3,512);
wave2=generate_wave_ellipsis_fm(20,120,pi/6,20,1,300,200,170,75,-pi/6,512);
patch1=generate_patch_rectangular(400,250,75,200,10,512);
patch2=draw_ellipsis(100,100,40,30,10,512);

wave=wave1+wave2+patch1+patch2;
