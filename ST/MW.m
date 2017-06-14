function [ tfr, as ] = MW( x, M, Ts, w0T, as_range, gamma_K, is_freq)
%[ tfr ] = MW( x, T, w0 )
%
% Compute the Morlet wavelet transform

if(~exist('gamma_K', 'var'))
 gamma_K = 10^(-5);
end

if ~exist('as_range', 'var')
 as_range = [0.05 1.5];
end

if ~exist('is_freq', 'var')
 is_freq = 0; 
end

as = a_axis(M,  as_range, is_freq); %logspace(log10(as_range(1)), log10(as_range(2)), M);  % scales

sqrt_pi = sqrt(pi);

x = x(:).';
N = length(x);

tfr = zeros(M, N);

for n = 1:N
 
 for m = 1:M
   K = round(sqrt(2 * log(1/gamma_K)) * (T/Ts) * as(m));
   

   k_min = min(n-1, round(K));
   k_max = min(N-n, round(K));
   
   k = (-k_min):k_max;

   tau = k * Ts;

   tfr(m,n) = Ts/(sqrt(abs(as(m)) * T * sqrt_pi)) * sum(x(n+k) .* exp(-tau.^2/(2*(T*as(m))^2)) .* exp(-1j * w0 * tau/as(m)));
 end
end

%tfr = transpose(tfr);

end
%% EOF