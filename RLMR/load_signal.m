% script load_signal
%
% generate a test signal s of length N, using variable signal = 1-6
%
% Tests signals:
%
% 1: Multicomponent signal made of 1 impulse, 1 sinsoid and 2 linear chirps
% 2: Multicomponent signal made of 3 impulses
% 3: Multicomponent signal made of 3 successive sinusoids at different frequencies
% 4: Multicomponent signal made of 3 simultaneous sinusoids at different frequencies
% 5: Elementary signal made of 1 pure sinusoid at 100Hz
% 6: Real audio speech signal
%
% Recommended requirements: TFTB (http://tftb.nongnu.org/index_fr.html)
%
% Author: D. Fourer (dominique@fourer.fr)
% Date: Sept 2015
% Ref: [D.Fourer, F. Auger and P.Flandrin. Recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. IEEE Proc. ICASSP 2016.]

if (signal == 1)
 Nchirp  =  420; loc_impulse1=10;loc_impulse2=50; val_impulse=15.0;  %%
 if ~exist('fmlin') && exist('signal1.mat', 'file')
     load('signal1.mat')
     warning('Please install TFTB (http://tftb.nongnu.org/index_fr.html), use default signal')
 else
 %% Oscillating frequency  chirp
%  t1 = (1:Nchirp)/Fs;
%  w1 = 2 * pi * 3;
%  a1 = 50;
%  p1 = pi;
%  freq_sin = 355;
%  f_inst = freq_sin + a1 *cos(w1 * t1 + p1);
%  phi_t = cumsum(f_inst) / Fs;
%  s0 = cos(2 * pi * phi_t );
%  
%   s = fmconst(Nchirp, 0.1) + fmlin(Nchirp,0.12,0.3) + s0';    % + fmlin(Nchirp,0.3,0.4);
 s = fmconst(Nchirp, 0.1) + fmlin(Nchirp,0.13,0.3) + fmsin(Nchirp,0.3,0.45,320,1,0.3,+1);
 N=500;
 s=[zeros(N-Nchirp,1);1.0*real(s)];
 s(loc_impulse1)=val_impulse;
 s(loc_impulse2)=val_impulse;
 end
elseif (signal == 12)  %% identical to 1 whithout the impulses
 Nchirp  =  500;
 if ~exist('fmlin') && exist('signal1.mat', 'file')
     load('signal1.mat')
     warning('Please install TFTB (http://tftb.nongnu.org/index_fr.html), use default signal')
 else
 %% Oscillating frequency  chirp
%  t1 = (1:Nchirp)/Fs;
%  w1 = 2 * pi * 3;
%  a1 = 50;
%  p1 = pi;
%  freq_sin = 355;
%  f_inst = freq_sin + a1 *cos(w1 * t1 + p1);
%  phi_t = cumsum(f_inst) / Fs;
%  s0 = cos(2 * pi * phi_t );
%  
%   s = fmconst(Nchirp, 0.1) + fmlin(Nchirp,0.12,0.3) + s0';    % + fmlin(Nchirp,0.3,0.4);
 s = fmconst(Nchirp, 0.1) + fmlin(Nchirp,0.13,0.3) + fmsin(Nchirp,0.3,0.45,320,1,0.3,+1);
 N=500;
 s = real(s);
 end
elseif (signal == 2)
 loc_impulse=15; val_impulse=10.0;
 N=1000; s=zeros(N,1); s(loc_impulse)=val_impulse; s(150)=val_impulse; s(300)=val_impulse;
 % m=117; nfreqs(m)*Fs, figure(6); plot(t,tfr(m,:),t,rtfr(m,:)); 

elseif (signal==3)
 N=1000; L1=350; L2=300; s=[real(fmconst(L1,0.1));real(fmconst(L2,0.2));real(fmconst(N-L1-L2,0.3))];
 % My_t=700; figure(6); plot(nfreqs,tfr(:,My_t),nfreqs,rtfr(:,My_t));

elseif (signal == 4)
 N=1000; s = real(fmconst(N, 0.1) + fmconst(N, 0.25) + fmconst(N, 0.4));
 % My_t=500; figure(6); plot(nfreqs*Fs,tfr(:,My_t),nfreqs*Fs,rtfr(:,My_t));

elseif (signal == 5)
  N=1000; s = real(fmconst(N, 0.1));
 % My_t=500; figure(6); plot(nfreqs*Fs,tfr(:,My_t),nfreqs*Fs,rtfr(:,My_t));

elseif (signal == 6)
  [s, Fs] = wavread('wav/voixr.wav');
  deb = 600;
  fin = deb + Fs-1;
  s = s(deb:fin);
  
  %% override user parameters 
  N = length(s);
  
  M = 4000;
  mi = 1;mf = 700; %550;
  k = 3; L = 30;
  mu_i= 0.8; mu_f=3;    %% mu range used for LM reassignment

end