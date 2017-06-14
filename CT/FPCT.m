function [Spec,f,t] = FPCT(Sig,SampFreq,Ratio,N,WinLen);

% Polynomial chirplet transform of freuqency domain
%
%	   Sig       : the signal to be analyzed
%    SampFreq    : sampling frequency
%      Ratio     : Amplitude of modulation component of the signal
%        N       : the number of frequency bins (default : length(Sig)).
%	WinLen       : the length of window used to locate the signal in time.
%
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

if(nargin < 3),
    error('At least 3 inputs are required!');
end

SigLen = length(Sig); % length of time series

if (nargin < 5),
    WinLen = N / 4;
end

if N ~= SigLen;
    N = SigLen;
end

RatioNum = length(Ratio);

%% compute frequency spectrum


Spec = fft(Sig)./length(Sig); %take the fft of our sin wave, y(t)
Freq=( -ceil((N-1)/2):N-1-ceil((N-1)/2) )/N*SampFreq;
df = Freq(end/2+1:end);
sigfft = Spec(1:end/2);

%% shift operator
shift = zeros(size(df));
for k = 1:RatioNum,
    shift = shift + Ratio(k) * df.^k;
end

%% rotation 
kernel = zeros(size(df));
for k = 1:RatioNum,
    kernel = kernel + Ratio(k)/(k+1) * df.^(k+1);
end

rSig = sigfft.* exp(j * 2 * pi * kernel');

WinLen = ceil(WinLen / 2) * 2;
f = linspace(-1,1,WinLen)';
WinFun = exp(log(0.005) * f.^2 );
WinFun = WinFun / norm(WinFun);

Lh = (WinLen - 1)/2; 

fft_len = length(df);
Spec = zeros(SigLen,fft_len) ;   % matrix
Rdt = zeros(SigLen,1);
conjSpec = Spec;

%% iteration
wait = waitbar(0,'Please wait...');
for iLoop = 1:fft_len,

    waitbar(iLoop/fft_len,wait);    
    tau = -min([round(fft_len/2)-1,Lh,iLoop-1]):min([round(fft_len/2)-1,Lh,fft_len-iLoop]);  % signal span
    temp = floor(iLoop + tau);  % temp is the signal covered by window
    temp1 = floor(Lh+1+tau);    % window span
    
    Rdt(1:length(temp)) = df(temp);         % time unit   
    wSig = rSig(temp);
    
    wSig = wSig.* conj(WinFun(temp1)); % Z(t)* complex conjugate of window?
    
    Spec(1:length(wSig),iLoop) = wSig;  % windowed analytic signal
    Spec(:,iLoop) = Spec(:,iLoop) .* exp(-j * 2.0 * pi * shift(iLoop) * Rdt);
    conjSpec(:,iLoop)  =  fliplr(Spec(:,iLoop));

end

iSpec  = [conjSpec(1:end-1,:);Spec];
iLen = length(iSpec);
Spec = iSpec(round(iLen/2):end,:);

Spec = ifft(Spec);
Spec = abs(Spec)/2/pi;

close(wait);

[SigLen,nLevel] = size(Spec);

f = [0:nLevel-1]/nLevel * SampFreq/2;  % frequency  in TF plane?
t = (0: SigLen-1)/SampFreq;      % time  in TF plane

[fmax fmin] = FreqRange(Sig);
fmax = fmax * SampFreq;
fmin = fmin * SampFreq;
Spec = Spec';
clf
%=====================================================%
% Plot the result                                     %
%=====================================================%
set(gcf,'Position',[20 100 350 300]);	    
set(gcf,'Color','w');					    

mesh(t,f,Spec); 
axis([min(t) max(t) fmin fmax]);
xlabel('Time / s');
ylabel('Freq / Hz');
Info = 'C = ';

for i = 1:RatioNum,
    Info = [Info,num2str(Ratio(i),4), '  '];
end

if RatioNum == 0,
    Info = 'C = 0';
end