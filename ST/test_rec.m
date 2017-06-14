clear all
close all

signal = 1; %5;  %% change value in [1,2,3,4]

Fs = 1000;          %% sampling frequency
Ts = 1/Fs;          %% sampling period

M_range = [100]; % 500; %[200 400];



w0T = 7 * pi; %7*pi;
gamma_K = 10^(-3);  %% 10^(-4)


%N=300;signal=5;                       %s = real(fmconst(N, 0.1));
load_signal
s = s-mean(s);              %% 0-mean signal


%%%%%%%%%%% 
% N = length(s);
% M_range = 100:10:500;
% q_results = zeros(1, length(M_range));
% for i = 1:length(M_range)
%     M = M_range(i);
% 
%     Fx = fft(s, M);
%     s_hat = N/M * ifft(Fx, N);
% 
% 
%     figure(1)
%     plot(real(s_hat))
%     hold on
%     plot(real(s), 'r-.')
%     hold off
%     q_results(i) = RQF(s, s_hat)
%     pause
% end
% 
% q_results
% pause
% Fx = fft(s, N);
% m = round(linspace(1, N, 400));
% m = m + 2;
% m = m(find(and(m>=1,m <=500)));
% s_hat = ifft(Fx(m), N);


for i = 1:length(M_range)
    
 M = M_range(i);
 fprintf(1, 'M=%d \n', M);
 
 %% version classique sans 0-padding
%  [tfr, ma] = tfrst(s, M, w0T, gamma_K);
%  s_hat = rectfrst(tfr, ma, M);
%  
%  fprintf(1, 'No 0-Padding RQF=%0.2f dB\n', RQF(s, s_hat'));
%  figure(1)
%  plot(s,'k-.')
%  hold on
%  plot(real(s_hat), 'r-.')
%  hold off
%


%  
%  %% version 1  0-padding
%  z_pad   = round(w0T * M * sqrt(2*log(1/gamma_K)) / pi);
%  s0p     = [zeros(z_pad,1) ; s ; zeros(z_pad, 1)];
%  i_s     = z_pad+1:(length(s0p)-z_pad);
%  %  %% Compute transfirm
%  [tfr0p] =  tfrst(s0p, M, w0T, gamma_K);
% 
%  %% reconstruction
%  stmp    = rectfrst(tfr0p);
%  s_hat   = stmp(i_s);
%  fprintf(1, '0-padding 1 RQF=%0.2f\n', RQF(s, s_hat'))
%   figure(2)
%   plot(s,'k-.')
%   hold on
%   plot(real(s_hat), 'r-.')
%  
%   
  
 if M >= N
  [tfr, ma] =  tfrstzp(s, M, w0T, gamma_K);
   %  %% reconstruction
   stmp   = rectfrst(tfr);  %, ma, M
   s_hat   = stmp(2:(end-1));
 else
   z_pad  = round(w0T * M * sqrt(2*log(1/gamma_K)) / pi);   %% maximal width
   s_tm = [zeros(1, z_pad)  s.'  zeros(1, z_pad)];
   i_s = (z_pad+1):(length(s_tm)-z_pad);
   [tfr, ma] =  tfrst(s_tm, M, w0T, gamma_K);
   stmp   = rectfrst(tfr, ma, M);
   s_hat  = stmp(i_s);
 end
 
 fprintf(1, '0-padding 2 RQF=%0.2f dB \n', RQF(s, s_hat.'));
  figure(3)
  plot(real(s),'k-.')
  hold on
  plot(real(s_hat), 'r-.')
  hold off

end





%%%%%%%% reconstruction
% mm = m_axis(M);
% I_nz = find(abs(mm) > eps);
% mm = mm(I_nz);
% tfr = tfr(I_nz, :);
% mm = transpose(mm);
% 
% %% 2 Compute Ch
% dK = 0.01;
% K  = 100;
% k =  ((w0T-K):dK:(w0T+K))+eps;
% C_h = dK * sum( 1./ k .* exp(-(k - w0T).^2 /2));
% 
% 
% R = zeros(1,N);
% 
% for im = 1:length(mm) %round(length(mm)/2)
%   for n = 1:N 
%      R(n) = R(n) + sign(mm(im)) * 1/C_h*  1./ mm(im) .* exp(1i * 2*pi*(n-1) * mm(im) / M) .* tfr(im, n); 
%   end
%   figure(1)
%   plot(real(R))
%   title(sprintf('m=%d [%d/%d], SNR=%.2f', mm(im), im, length(mm), SNR(s, R')));
%   pause
%   %R(n) = sum(1./abs(mm) .* exp(1i * 2*pi*(n-1) * mm / M) .* tfr(:, n));
% end
% 
% % for n = 1:N 
% %   R(n) = sum(1./abs(mm) .* exp(1i * 2*pi*(n-1) * mm / M) .* tfr(:, n));
% % end
% 
% 
% %s_hat = 1 / C_h * R ;