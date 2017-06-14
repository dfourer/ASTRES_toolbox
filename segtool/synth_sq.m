function [s] = synth_sq(Tx,dt,mywav,nv)
%SYNTH_SQ synthesis

if nargin<4
    mywav = 'cmor2-1';
    nv=32;
end

N = size(Tx,2);
t = dt*(0:N-1);
noct = log2(N)-1;
na = noct*nv;
as = (2^(1/nv) .^ (1:1:na));
assert(na==size(Tx,1));


%% Calcul de CPsi
PsiMor=@(x,Fc,Fb) 0.5*exp(-pi^2*Fb*(x-Fc).^2) ./x;
if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    beta = str2num(mywav(v1:v2-2));
    gamma = str2num(mywav(v2:end));
    PsiMor=@(x) 0.5*gmor(beta,gamma,x)./x;
    CPsi = quadgk(@(x)PsiMor(x),eps,(beta/gamma)^(1/gamma)+2*sqrt(gamma*beta));  
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    Fb = str2num(mywav(v1:v2-2));
    Fc = str2num(mywav(v2:end));
    PsiMor=@(x) 0.5*cmor(Fb,Fc,x)./x;
    CPsi = quadgk(@(x)PsiMor(x),0,inf);  
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9]*-[0-9]');
    mu = str2num(mywav(v1:v2-2));
    sigma = str2num(mywav(v2:end));
    PsiMor=@(x) 0.5*bump(mu,sigma,x)./x;
    CPsi = quadgk(@(x)PsiMor(x),0,inf);  
end 


%% Somme pondérée de Tx;
s = real(sum(Tx) / CPsi);
s = s/nv*log(2);


end

