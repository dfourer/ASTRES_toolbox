function [fmax, fmin] = FreqRange(Sig);
%======================================================
% To decide the frequency range used for calculating     
% the PTFR                         
% 
%     Sig           any one dimension sig         
%    fmax           the maximal frequency  
%    fmin           the minimum frequency
%
% Z K Peng
% Pim, Tsinghua University 
% 02, Sept, 2001
% All rights reserved 
%======================================================

if(length(Sig) < 32),
    error('The signal is too short!');
end
sptr = fft(Sig);                                      %sptr is the FFT spectrum of the signal 
sptrLen = length(sptr);
sptr = sptr(1:round(sptrLen/2)+1);
f = linspace(0,0.5,round(sptrLen/2)+1);

iLen = length(sptr);
TotalEngery = sum(abs(sptr));
MaxSptr = max(abs(sptr)); 

Temp = 0;
maxFreq = 0;
for j = 1:iLen
    Temp = Temp + abs(sptr(j));
    if(Temp > 0.90 * TotalEngery)
        maxFreq = j;
        break;
    end
end

Temp = find(abs(sptr) > MaxSptr/50);
fmin = f(Temp(1));

Temp = max(Temp);
if(Temp > maxFreq), 
    maxFreq = Temp;
end

Temp = f(maxFreq);
fmax = ceil(Temp/0.05) * 0.05;

fmin = floor(fmin/0.05) * 0.05;
if(fmin == 0.0),
    fmin = 0.001;
end