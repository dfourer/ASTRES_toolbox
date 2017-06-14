function FPCT_example(),
% example to demonstrate the use of FPCT

%   Written by Yang Yang, March 2011.  
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
clear all
close all

SampFreq = 100;
t = 0:1/SampFreq:15;
Sig = exp(j*(-1*2*pi*t - 10*pi*t.^2 + 0.2*2*pi*t.^3 ));
Sig = [fliplr(Sig(2:end)) Sig];

iSig = ifft(Sig);
SigLen = length(iSig);
iSig = iSig(1:round(SigLen/2));

figure(1)
plot(t,abs(iSig))
set(gcf,'Position',[20 100 350 300]);	    
set(gcf,'Color','w');

figure(2)
[Spec,f] = FPCT(iSig',SampFreq,0,1024,256); %SFFT

%% 
[v, I] = max(Spec,[],2);
[p, z] = polylsqr(f,t(I),2);
figure(3)
FPCT(iSig',SampFreq,z(2:end),1024,256); %STFT
z