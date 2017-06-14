function [Wx Tx fs as] = sqt(s,dt,gamma,mywav,nv,doplot)
% function sqt : computes the synchrosqueezing transform of signal s.
% Implements algorithm of Brevdo, Daubechies, Wu et al.
%
% Inputs
%   s : input signal, power of 2
%   dt : sample period
%   gamma : threshold
%   mywav : string coding the mother wavelet (see mycwt)
%   nv : number of coefficient per octave
%   doplot : 0 or 1 for extern plotting
%
% Outputs : 
%   Wx : the wavelet transform
%   Tx : the synchrosqueezed transform
%   fs : frequency vector
%   as : scales vector


% Parameters
if nargin<5
    nv = 32; % nbre coef    s par octave
    gamma = 0.01; % Seuil pour wavelet thresholding
    mywav = 'cmor2-1';
end

if nargin<6
    doplot=0;
end


% Calculate the wavelet transform - padded via symmetrization
s = s(:);
N = length(s);
t = dt*(0:N-1);

noct = log2(N)-1;
assert(noct > 0 && mod(noct,1) == 0);
assert(nv>0 && mod(nv,1)==0); % integer
assert(dt>0);
assert(~any(isnan(s)));
na = noct*nv;
as = (2^(1/nv) .^ (1:1:na));

% Compute continous wavelet transform
Wx = mycwt(s,as,mywav,dt);

% Candidate instantaneous frequency
Wxx = unwrap(angle(Wx),[],2);
Wxr = [Wxx(:, end-1:end) Wxx Wxx(:, 1:2)];
w2 = -Wxr(:, 5:end);
w2 = w2 + 8*Wxr(:, 4:end-1);
w2 = w2 - 8*Wxr(:, 2:end-3);
w2 = w2 + Wxr(:, 1:end-4);
w2 = w2 / (12*dt*2*pi);
w2(abs(Wx)<gamma)=NaN;
w2(abs(w2)>N/2) = NaN;
w2(w2<0) = NaN;
w = w2;


% Synchrosqueezing transform
fs = as;
Tx = zeros(size(Wx));
for b=1:N
    for ai=1:length(as)
        if ~isnan(w(ai,b)) && abs(Wx(ai,b))>gamma
        	k = floor(1+nv*log2(w(ai,b)));
            if k>0 && k<=na
            	Tx(k,b) = Tx(k, b) + Wx(ai, b);
            end
        end
    end % for ai...
end % for b


% Plotting
if doplot
    sqplot(Tx,dt,fs);
end

end

