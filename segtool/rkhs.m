function Kab = rkhs(mywav,dt,as,N,a,b)
% Kab = computes the reproducing kernel associated to the mother wavelet
% mywav. test function, not used in segtool.

na = length(as);
xi = dt*[0:N/2 -N/2+1:-1];

Kab = zeros(na,N);

% Selecting mother wavelet
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

%figure();plot((real(ifft((filt(as(100)))))));return

psia = filt(as(a));
for ap=1:na
    Kab(ap,:) = as(ap)*fftshift(ifft(psia .* conj(filt(as(ap)))));
end

Kab = circshift(Kab,[0 -N/2+b-1]);
