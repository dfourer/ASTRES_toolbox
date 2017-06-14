% select the test signal using variable signal = 1-6
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
% Author: D. Fourer (dominique@fourer.fr)
% Date: July 2015
% Ref: [D.Fourer and F. Auger. Real-time recursive versions of the Levenberg-Marquardt
% reassigned spectrogram and of the synchrosqueezed STFT. Proc. ICASSP 2016.]
if ~exist('complex_signal', 'var')
 complex_signal = 0;
end

if (signal == 1)
 Nchirp  =  440; loc_impulse1=15;loc_impulse2=40; val_impulse=10.0;  %%
 
 %% Oscillating frequency  chirp
 t1 = (1:Nchirp)/Fs;
 w1 = 2 * pi * 3;
 a1 = 50;
 p1 = pi;
 freq_sin = 355;
 f_inst = freq_sin + a1 *cos(w1 * t1 + p1);
 phi_t = cumsum(f_inst) / Fs;
 s0 = cos(2 * pi * phi_t );
 
 s = fmconst(Nchirp, 0.1) + fmlin(Nchirp,0.12,0.3) + s0';    % + fmlin(Nchirp,0.3,0.4);
 N=500;
 s=[zeros(N-Nchirp,1);1.0*s];
 s(loc_impulse1)=val_impulse;
 s(loc_impulse2)=val_impulse;
 
 if complex_signal == 0   %% otherwise s is complex with different real / imaginary parts
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
  N=500; %s = real(fmconst(N, 0.1));
  
  s = fmconst(N, 0.1);
  %s = 2 * exp(1i * 2*pi * 0.1 * (0:(N-1)));
  

  %% impulse
  %N=300; 
  %s = zeros(N,1);
  %s(100)  = 20;
  
  %% const sine
  %s = fmconst(N, 0.1);
  
  %chirp
  %s = fmlin(N,0.12,0.3);
  
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