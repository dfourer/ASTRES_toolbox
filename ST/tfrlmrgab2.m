function [tfr,rtfr,hat] = tfrlmrgab(x,t,N,mu,Nh,trace,K)
%TFRRGAB Reassigned Gabor spectrogram time-frequency distribution.
%	[TFR,RTFR,HAT] = TFRLMRGAB(X,T,N,NH,TRACE,K) 
%	computes the Gabor spectrogram and its reassigned version.
%	This particular window (a Gaussian window) allows a 20 % faster
%	algorithm than the TFRRSP function.
%                  
%	X     : analysed signal
%	T     : the time instant(s)           (default : 1:length(X))
%	N     : number of frequency bins      (default : length(X))
%	NH    : length of the gaussian window (default : N/4))
%	TRACE : if nonzero, the progression of the algorithm is shown
%                                         (default : 0).
%	K     : value at both extremities     (default 0.001)
%	TFR,  : time-frequency representation and its reassigned
%	RTFR    version. When called without output arguments, 
%	        TFRRGAB runs TFRQVIEW.
%	HAT   : Complex matrix of the reassignment vectors.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); tfrrgab(sig,1:128,128,19,1);
%
%	See also all the time-frequency representations listed in
%	 the file CONTENTS (TFR*)

%	F. Auger, May-July 1994, July 1995. 
%       Copyright (c) 1996 by CNRS(France). 
% 
%	------------------- CONFIDENTIAL PROGRAM --------------------
% 	This program can not be used without the authorization of its 
% 	author(s). For any comment or bug report, please send e-mail to 
% 	f.auger@ieee.org

if (nargin == 0),
 error('At least 1 parameter required');
end;
[xrow,xcol] = size(x);
if (nargin <= 2),
 N=xrow;
end;

hlength=floor(N/4);
hlength=hlength+1-rem(hlength,2);

if (nargin == 1),
 t=1:xrow; 
end;

if (nargin <= 4),
 Nh=hlength; trace=0; K=0.001;
elseif (nargin == 5),
 trace = 0; K=0.001;
elseif (nargin == 6),
 K= 0.001;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N & nargin==6),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

if (rem(Nh,2)==0), 
 error('Nh must be odd'); 
elseif length(Nh)~=1,
 error('Nh must be a scalar');
end;

Nh2=Nh-2;
TFTBcontinue=1;
while TFTBcontinue,
 Nh2=Nh2+2;
 h=tftb_window(Nh2,'gauss',K^((Nh2-1)^2 /(Nh-1)^2)); 
 TFTBcontinue=(h(Nh2)*(Nh2-1)>2*K);
end;

K=K^((Nh2-1)^2 /(Nh-1)^2); Nh=Nh2; Lh=(Nh-1)/2; 
h=h; 
Th=h.*[-Lh:Lh]'; T2h=Th.*[-Lh:Lh]';

if (tcol==1),
 Dt=1; 
else
 Deltat=t(2:tcol)-t(1:tcol-1); 
 Mini=min(Deltat); Maxi=max(Deltat);
 if (abs(Mini-Maxi) > eps), %Mini~=Maxi
  error('The time instants must be regularly sampled.');
 else
  Dt=Mini;
 end;
 clear Deltat Mini Maxi;
end;

tfr= zeros(N,tcol); tf2= zeros(N,tcol); tf3= zeros(N,tcol); tf4= zeros(N,tcol);
if trace, disp('Gabor spectrogram'); end;

for icol=1:tcol,
 if trace, DISPROG(icol,tcol,10); end;
 ti= t(icol); 
 tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1;
 norm_h=norm(h(Lh+1+tau));
 tfr(indices,icol)=x(ti+tau).*conj(  h(Lh+1+tau))/norm_h;
 tf2(indices,icol)=x(ti+tau).*conj( Th(Lh+1+tau))/norm_h;
 tf4(indices,icol)=x(ti+tau).*conj(T2h(Lh+1+tau))/norm_h; 
end ;
tfr=fft(tfr); tf2=fft(tf2); tf4=fft(tf4); 

alpha=-log(K)/Lh^2;

avoid_warn=find(tfr~=0.0);
tf3(avoid_warn)=imag(-2*alpha*tf2(avoid_warn)./tfr(avoid_warn));
tf2(avoid_warn)=tf2(avoid_warn)./tfr(avoid_warn);
tf4(avoid_warn)=tf4(avoid_warn)./tfr(avoid_warn);
tfr=abs(tfr).^2;
if trace, fprintf ('\nLevenberg-Marquardt reassignment: \n'); end;

rtfr= zeros(N,tcol); 
Ex=mean(abs(x(min(t):max(t))).^2); Threshold=1.0e-6*Ex;
for icol=1:tcol,
 if trace, disprog(icol,tcol,10); end;
 for jcol=1:N,
  if abs(tfr(jcol,icol))>Threshold,
   icolhat= icol + round(real(tf2(jcol,icol)/Dt)); icolhat=min(max(icolhat,1),tcol);
   jcolhat= jcol - round(tf3(jcol,icol)*N/(2.0*pi)); jcolhat=rem(rem(jcolhat-1,N)+N,N)+1;
   
   d2Phi_dw2= -imag(tf4(jcol,icol)-tf2(jcol,icol)^2);
   d2Phi_dt2= -4*alpha^2*d2Phi_dw2;
   d2Phi_dtdw=-2*alpha*real(tf4(jcol,icol)-tf2(jcol,icol)^2);
   a11=1+mu+d2Phi_dtdw; a12=d2Phi_dw2; a21=-d2Phi_dt2; a22=mu-d2Phi_dtdw;
   Det=a11*a22-a12*a21;
   r1=-real(tf2(jcol,icol)); r2=tf3(jcol,icol);
   icoltilde=icol-round((a22*r1-a12*r2)/Det/Dt);icoltilde=min(max(icoltilde,1),tcol);
   jcoltilde=jcol+round((a21*r1-a11*r2)*N/(2.0*pi*Det));jcoltilde=rem(rem(jcoltilde-1,N)+N,N)+1;
   %while (jcolhat<1),jcolhat=jcolhat+N; end;
   %while (jcolhat>N),jcolhat=jcolhat-N; end;
   rtfr(jcoltilde,icoltilde)=rtfr(jcoltilde,icoltilde) + tfr(jcol,icol) ;
   tf2(jcol,icol)=jcoltilde + 1j * icoltilde;
   %fprintf('icol=%d ; jcol=%d ; d2Phi_dw2=%f ; d2Phi_dt2=%f; d2Phi_dtdw=%f; Det=%f; r1=%f; r2=%f; icoltilde=%d, jcoltilde=%d\n', ...
   %        icol, jcol, d2Phi_dw2, d2Phi_dt2, d2Phi_dtdw, Det, real(r1), r2, icoltilde, jcoltilde);
  else
   tf2(jcol,icol)=inf*(1+j);
   rtfr(jcol,icol)=rtfr(jcol,icol) + tfr(jcol,icol) ;
  end;
 end;
end;

if trace, fprintf('\n'); end;
clear tf3;
if (nargout==0),
 TFTBcontinue=1;
 while (TFTBcontinue==1),
  choice=menu ('Choose the representation:',...
               'stop',...
               'Gabor spectrogram',...
               'LM reassigned Gabor spectrogram');
  if (choice==1), TFTBcontinue=0;
  elseif (choice==2), 
   Q=round(tcol*N/xrow);
   tfrqview(tfr,x,t,'tfrgabor',tcol,Q,h);
  elseif (choice==3),
   tfrqview(rtfr,x,t,'type1',Nh);
  end;
 end;
elseif (nargout>2),
 hat=tf2;
end;

