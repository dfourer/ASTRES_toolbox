function [Wx,  as, stfr, f, lost] = sMW( x, M, Ts, T, w0, as_range, gamma_K, is_freq)
% [Wx,  as, stfr, f, lost] = my_sMW( x, M, Ts, T, w0, as_range, gamma_K )
%
% Compute the synchrosqueezed Morlet wavelet transform
%
% INPUT:
% x        : the signal to process
% M        : number of frequency bins to process (default: length(x))
% Ts       : the sampling period
% T        : T parameter of the Morlet wavelet transform
% w0       : central frequency in radian of the Morlet wavelet transform
% as_range : scale range to compute (defautl is [w0/(2*pi*50) w0/(pi/Ts)] and corresponds to [50Hz -- Fs/2Hz])
% gamma_K  : threshold applied on window values (default: 10^(-4))
%
% OUTPUT:
% Wx     : discrete Morlet wavelet transform
% as     : scales associated to y coordinates
% stfr   : synchrosqueezed Morlet wavelet transform
% f      : frequencies associated to y coordinates
% lost   : amount of lost energy during reassignment
%
% Author: D.Fourer (dominique@fourer.fr)
% Date: 28-08-2015
% Ref: [D. Fourer, F. Auger and J. Hu, Reassignment and Synchrosqueezing of the Stockwell Transform. IEEE. Trans. Signal Processing. 2015]
% Ref: [I. Daubechies, J. Lu and H-T Wu, Synchrosqueezed wavelet transforms: An empirical mode decomposition-like tool,
%       Applied and Computational Harmonic Analysis. vol 30(2), pages 243-261. 2011]

if(~exist('gamma_K', 'var'))
 gamma_K = 10^(-4);
end

if ~exist('as_range', 'var')
  as_range = [w0/(2*pi*50) w0/(pi/Ts)];    %% [50Hz -- Fs/2Hz]
end

if ~exist('is_freq', 'var')
 is_freq = 0; 
end

x = x(:).';   %%  convert to column vector
N = length(x);
lost = 0;
%% compute scales
as = a_axis(M,  as_range, is_freq); %as = logspace(log10(as_range(1)), log10(as_range(2)), M);

%% adjust frequencies
f = w0 ./ (2*pi*as);                         %% convert scales to frequencies in Hz
I = find(and(~isnan(f),f <= 1/(2*Ts)));      %% choose correct frequencies

as = as(I);
f = w0 ./ (2*pi*as);

stfr = zeros(M, N);
Wx    = zeros(M, N);
Wx_d  = zeros(M, N);

sqrt_pi = sqrt(pi);

for n = 1:N
 for m = 1:M
   K = round(sqrt(2 * log(1/gamma_K)) * (T/Ts) * as(m));
   k_min = min(n-1, round(K));
   k_max = min(N-n, round(K));
   k = (-k_min):k_max;
   tau = k * Ts; 
   A_m = 1/(sqrt(abs(as(m)) * T * sqrt_pi));
   
   if m > 1
     ds = as(m)-as(m-1);
   else
     ds = as(2)-as(1);  %% assume the first delta_s is equal to the second one
   end
   
   Wx(m,n)   = Ts * A_m * sum(x(n+k) .* conj(psi(tau/as(m), T, w0)));   %% MW_x(t,s)
   tmp_tfr = abs(Wx(m,n))^2;
   
   if tmp_tfr > eps       
    Wx_d(m,n) = -Ts * A_m * (T*as(m))^(-2) * sum(x(n+k) .* tau .* conj(psi(tau/as(m), T, w0))) - 1j*w0/as(m) * Wx(m,n);   %% MW_x^Dh(t,s)
    
    f_hat     = -imag(Wx_d(m,n) / Wx(m,n)) / (2*pi);

    [~, m_hat] = min(abs(f-f_hat));   %% match frequencies
    
    %% out of bounds
    m_out_of_bounds = m_hat < 1 || m_hat > M;
      
    if m_out_of_bounds
     lost = lost + tmp_tfr;
      continue;
    end
   
    stfr(m_hat, n) = stfr(m_hat, n) + Wx(m,n) * as(m)^(-3/2) * abs(ds);
   end
 end %% M
end %% N

% Wx   = transpose(Wx);
% stfr = transpose(stfr);
end

function v = psi(t, T, w0)
 v = exp(-t.^2/(2*T^2)) .* exp(1j * w0 * t);
end

%% EOF
