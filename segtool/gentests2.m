function [t dt x2cos x2chirpl x2chirpexp mchirp] = gentests2()
% gentests2 : generates test signals for SQ numerical tests


N = 1024;
dt = 1/1024;
t = dt*(0:N-1);

x2cos = cos(2*100*pi*t)+1.5*sin(2*129*pi*t);

mchirp = cos(2*77*pi*(t+0.2).^2);

x2chirpl = cos(2*250*pi*t.^2)+1.5*sin(2*329*pi*t.^2);

x2chirpexp = cos(2*20*pi*exp(1.5*t))+1.5*sin(2*29*pi*exp(1.5*t));

end