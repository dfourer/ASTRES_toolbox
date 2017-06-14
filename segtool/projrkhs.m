function Wtilde = projrkhs(Wx,mywav,dt,as,N,nv)
% projrkhs = computes the projection onto the RKHS WL2 of the truncated
% wavelet transform Wx.
%
% Inputs
%   Wx : any truncated wavelet transform
%   mywav : string, name of the mother wavelet
%   dt : sample period
%   as : scale vector
%   N : number of time samples
%   nv : number of coefficients per octave
% Output
%   Wtilde : projection of Wx onto WL2

na = length(as);
xi = 2*dt*[0:N/2 -N/2+1:-1];

Kab = zeros(na,N);
K = Kab;
Wtilde = zeros(na,N);

% Selecting mother wavelet
if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    beta = str2num(mywav(v1:v2-2));
    gamma = str2num(mywav(v2:end));
    filt = @(a) gmor(beta,gamma,a*xi);
    PsiMor=@(x) abs(gmor(beta,gamma,x)).^2 ./x;
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    Fb = str2num(mywav(v1:v2-2));
    Fc= str2num(mywav(v2:end));
    filt = @(a) cmor(Fb,Fc,a*xi);
    PsiMor=@(x) abs(cmor(Fb,Fc,x)).^2 ./x;
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    mu = str2num(mywav(v1:v2-2));
    sigma= str2num(mywav(v2:end));
    filt = @(a) bump(mu,sigma,a*xi);
    PsiMor=@(x) abs(bump(mu,sigma,x)).^2 ./x;
end

Cpsi = quadgk(PsiMor,0,inf);  
Delta = 2*floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps)+2;

% Classical implementation
for a=1:na
    K = zeros(na,N);
    psia = filt(as(a));
    %for ap=1:na
    for ap = max(a-Delta,1):min(na,a+Delta)
        phi = ifft(psia .* conj(filt(as(ap))));
        K(ap,:) = phi';
    end
    %K = K .* (as(:).^(0) * ones(1,N));
    for b=1:N
        Kab = circshift(K,[0 b-1]);
        if b==-500
            imagesc(abs(Kab)); pause;
        end
        Wtilde(a,b) = sum((Wx(:) .* Kab(:)));
    end
end

Wtilde = Wtilde / Cpsi  / nv / log(2) /2;

figure();imagesc(log(1+abs(Wtilde)));
figure();imagesc(log(1+abs(Wx)));


