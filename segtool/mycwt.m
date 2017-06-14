function Wx = mycwt(s,as,mywav,dt)
% mycwt : computes the continuous wavelet transform of a signal
%
% Inputs :
%   s : signal, power of 2
%   as : vector of scales
%   mywav : string name of the wavelet
%   dt : sample period
% Output
%   Wx : matrix of wavelet coefficients. Wx(a,b) : coefficient at scale
%   as(a) and time b*dt.

s = s(:); 
n = length(s);
na = length(as);


% symmetric padding
[N x n1] = mypad(s);

% Frequency vector : fftshift
xi = dt*[0:N/2 -N/2+1:-1];

Wx = zeros(na, N);
x = x(:).';
xh = fft(x);

% Filter definition
if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    beta = str2num(mywav(v1:v2-2));
    gamma = str2num(mywav(v2:end));
    filt = @(a) gmor(beta,gamma,a*xi);
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    Fb = str2num(mywav(v1:v2-2));
    Fc= str2num(mywav(v2:end));
    filt = @(a) cmor(Fb,Fc,a*xi);
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    mu = str2num(mywav(v1:v2-2));
    sigma= str2num(mywav(v2:end));
    filt = @(a) bump(mu,sigma,a*xi);
end


% for each octave
for ai = 1:na
	a = as(ai);
    psih = conj(filt(a));
    Wx(ai, :) = ifft(psih .* xh);
end

Wx = Wx(:, n1+1:n1+n);

end 
