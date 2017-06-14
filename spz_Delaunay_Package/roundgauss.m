function [w, L]=roundgauss(n,k);
% [W, L]=ROUNDGAUSS(N,PREC) compute the Gaussian
% window that is circular in the time-frequency plane.
%
%	N : number of frequency bins,
%	PREC : value of window at boundaries (default: PREC=1e-3)
%	L :  window length
%	W :  Gauss window.
%
% Example :
% w=roundgauss(512); plot(w); 
%
% Proof :
% by def. : w(t)=exp(-Pi*(t/T)^2)
% circular are obtained for T=sqrt(N)
% we set w(L/2)=PREC => L=2*sqrt(-N/pi log(PREC))
%

% E. Chassande-Mottin, February 1997

if nargin==1,
 k=1e-3;
end;

L=sqrt(n);
l=floor(sqrt(-n*log(k)/pi))+1;
w=amgauss(2*l+1,l+1,L);
w=w/norm(w);