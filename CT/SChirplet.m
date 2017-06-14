function [Spec,f] = SChirplet(Sig,SampFreq,Ratio,shapoint,N,WinLen);
% Cubic spline kerneled chirplet Transform
%
%	   Sig       : the signal to be analyzed
%    SampFreq    : sampling frequency
%      Ratio     : coeff. of each piece of low order polynomial  
%     shapoint   : Shape point of spline kernel
%        N       : the number of frequency bins (default : length(Sig)).
%	WinLen       : the length of window used to locate the signal in time.

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


if(nargin < 4),
    error('At least 4 inputs are required!');
end

SigLen = length(Sig);

if (nargin < 6),
    WinLen = SigLen / 4;
end

if (nargin < 5),
    N = SigLen;
end

if length(Ratio)==1
    Ratio = zeros(1,4);
end
  

    
sig = Sig;

RatioNum = size(Ratio,1);

dt = (0:(SigLen-1))';
dt = dt/ SampFreq;
if length(shapoint)==1
    shapoint = [1 dt(end)];
end

sf = zeros(size(dt));
kernel = zeros(size(dt));  

toMod = 0;
for k = 1:RatioNum,
        xi = shapoint(k);
        xii = shapoint(k+1);
        loc = find(dt>=xi & dt<=xii);
        ti = dt(loc)-xi;
        a = Ratio(k,1);
        b = Ratio(k,2);
        c = Ratio(k,3);
        d = Ratio(k,4);   
        kernel(loc) = (a*ti.^4/4 + b*ti.^3/3 + c*ti.^2/2 + d *ti)*2*pi;  %rotate
        kernel(loc) = kernel(loc)+toMod;  
        toMod = kernel(loc(end));            
        sf(loc) = (a*ti.^3 + b*ti.^2 + c*ti + d); % shift
end
% plot(dt,kernel)
% pause

rSig = hilbert(real(Sig));  %Z(t)
Sig = rSig .* exp(-j*kernel);

WinLen = ceil(WinLen / 2) * 2;
t = linspace(-1,1,WinLen)';
WinFun = exp(log(0.005) * t.^2 );
WinFun = WinFun / norm(WinFun);
Lh = (WinLen - 1)/2; 

Spec = zeros(N,SigLen) ;   % matrix

wait = waitbar(0,'Please wait...');
for iLoop = 1:SigLen,

    waitbar(iLoop/SigLen,wait);
    
    tau = -min([round(N/2)-1,Lh,iLoop-1]):min([round(N/2)-1,Lh,SigLen-iLoop]);  % signal span
    temp = floor(iLoop + tau);
   
    rSig = Sig .* exp(j*2*pi*sf(iLoop)*dt); % shift: IF
    rSig = rSig(temp);
    
    temp1 = floor(Lh+1+tau);    % window span
    rSig = rSig .* conj(WinFun(temp1)); % Z(t)* complex conjugate of window?
    Spec(1:length(rSig),iLoop) = rSig;  % windowed analytic signal
end;

Spec = fft(Spec); 
Spec = abs(Spec);

close(wait);

Spec = Spec(1:round(end/2),:);
[nLevel, SigLen] = size(Spec);

f = [0:nLevel-1]/nLevel * SampFreq/2;  % frequency  in TF plane?
t = (0: SigLen-1)/SampFreq;      % time  in TF plane

[fmax fmin] = FreqRange(sig);
fmax = fmax * SampFreq;
fmin = fmin * SampFreq;

clf
%=====================================================%
% Plot the result                                     %
%=====================================================%
set(gcf,'Position',[20 100 350 300]);	    
set(gcf,'Color','w');					    

mesh(t,f,Spec);  
colormap jet;
shading interp;
axis([min(t) max(t) fmin fmax]);
% colorbar;
ylabel('Freq / Hz');
xlabel('Time / Sec')
Info = 'C = ';

for i = 1:RatioNum,
    Info = [Info,num2str(Ratio(i),4), '  '];
end

if RatioNum == 0,
    Info = 'C = 0';
end

%title(['Nonlinear Chirplet[',Info,']']);
