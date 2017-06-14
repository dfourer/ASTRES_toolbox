function SChirplet_example()

% example to demonstrate the use of SChirplet

%   Written by Yang Yang , March 2011.
%   email: z.peng@sjtu.edu.cn
%	Copyright (c) belongs to the authors of the papers: 
%)	Peng Z.K , Meng G., Lang Z.Q.,Chu F.L, Zhang W.M., Yang Y., Polynomial Chirplet Transform with Application to Instantaneous Frequency Estimation,
%   IEEE Transactions on Measurement and Instrumentation 60(2011) 3222-3229
%	Yang Y., Peng Z.K., Dong X.J., Zhang W.M., Meng G.,General parameterized time-frequency transform, 
%   IEEE Transactions on Signal Processing, 62(2014) 2751-2764
%   Yang Y., Peng Z.K., Zhang W.M., Meng G., Spline-kernelled chirplet transform for the analysis of signals with time-varying frequency and its application,
%   IEEE Transactions on Industrial Electronics, 59(2012) 1612-1621
%   Yang Y., Peng Z.K, Meng G., Zhang W.M., Characterize highly oscillating frequency modulation using generalized Warblet transform,
%   Mechanical Systems and Signal Processing, 26 (2012) 128-140
%   Yang Y., Peng Z.K., Zhang W.M., Meng G., Frequency-varying group delay estimation using frequency domain polynomial chirplet transform, 
%   Mechanical Systems and Signal Processing, 46(2014) 146-162
%   The citation about the papers must be included in all publications or
%   thesises as long as this program is used by anyone.


clc
clear
close all
SampFreq = 100;
t = 0 : 1/SampFreq : 10;
Sig = sin(2*pi*(10*t+5*atan((t - 5).^2)));

figure(1)
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w'); 
plot(t,Sig);
axis tight
xlabel('Time/Sec');

figure(2)
set(gcf,'Position',[20 100 320 250]);	    
set(gcf,'Color','w');	
Spec = 2*fft(Sig)/length(t);              % FFT
Freq = linspace(0,SampFreq, length(t));
plot(Freq(1:end/2),Spec(1:end/2));
axis tight
xlabel('Frequency / Hz');

N = 512;
WinLen =128;
figure(3)
[Spec,f] = SChirplet(Sig',SampFreq,0,1,N,WinLen);

p = 25; % piece number of cubic spline
[v l] = max(Spec,[],1);
sp = splinefit(t,f(l),p);
pp = ppval(sp,t);

figure(4)
Ratio = sp.coefs;
shapepoint = sp.breaks;
WinLen =1024;
[Spec,f] = SChirplet(Sig',SampFreq,Ratio,shapepoint,N,WinLen);