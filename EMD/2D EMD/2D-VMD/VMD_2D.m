function [u, u_hat, omega] = VMD_2D(signal, alpha, tau, K, DC, init, tol)
% 2D Variational Mode Decomposition
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% {konstantin,zosso}@math.ucla.edu
% http://www.math.ucla.edu/~{konstantin,zosso}
% Initial release 2014-03-17 (c) 2014
%
% Input and Parameters:
% ---------------------
% signal     - the space domain signal (2D) to be decomposed
% alpha      - the balancing parameter for data fidelity constraint
% tau        - time-step of dual ascent ( pick 0 for noise-slack )
% K          - the number of modes to be recovered
% DC         - true, if the first mode is put and kept at DC (0-freq)
% init       - 0 = all omegas start at 0
%              1 = all omegas start initialized randomly
% tol        - tolerance of convergence criterion; typically around 1e-7
%
% When using this code, please do cite our papers:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing, 62(3):531-544, 2014. DOI:10.1109/TSP.2013.2288675
%
% K. Dragomiretskiy, D. Zosso, Two-Dimensional Variational Mode
% Decomposition, IEEE Int. Conf. Image Proc. (submitted). Preprint
% available here: ftp://ftp.math.ucla.edu/pub/camreport/cam14-16.pdf
%


% Resolution of image
[Hy,Hx] = size(signal);
[X,Y] = meshgrid((1:Hx)/Hx, (1:Hy)/Hy);


% Spectral Domain discretization
fx = 1/Hx;
fy = 1/Hy;
freqs_1 = X - 0.5 - fx;
freqs_2 = Y - 0.5 - fy;

% N is the maximum number of iterations
N=3000;

% For future generalizations: alpha might be individual for each mode
Alpha = alpha*ones(K,1);

% Construct f and f_hat
f_hat = fftshift(fft2(signal));

% Storage matrices for (Fourier) modes. All iterations are not recorded.
u_hat = zeros(Hy,Hx,K);
u_hat_old = u_hat;
sum_uk = 0;

% Storage matrices for (Fourier) Lagrange multiplier.
mu_hat = zeros(Hy,Hx);

% N iterations at most, 2 spatial coordinates, K clusters
omega = zeros(N, 2, K);

% Initialization of omega_k
switch init
    case 0
        % spread omegas radially
        
        % if DC, keep first mode at 0,0
        if DC
            maxK = K-1;
        else
            maxK = K;
        end
        for k = DC + (1:maxK) 
            omega(1,1,k) = 0.25*cos(pi*(k-1)/maxK);
            omega(1,2,k) = 0.25*sin(pi*(k-1)/maxK);
        end
        
        % Case 1: random on half-plane
    case 1
        for k=1:K
            omega(1,1,k) = rand()-1/2;
            omega(1,2,k) = rand()/2;
        end
        
        % DC component (if expected)
        if DC == 1
            omega(1,1,1) = 0;
            omega(1,2,1) = 0;
        end
end

%% Main loop for iterative updates

% Stopping criteria tolerances
uDiff=tol+eps;
omegaDiff = tol+eps;

% first run
n = 1;

% run until convergence or max number of iterations
while( ( uDiff > tol || omegaDiff > tol ) && n < N )
    
    % first things first
    k = 1;
    
    % compute the halfplane mask for the 2D "analytic signal"
    HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
    
    % update first mode accumulator
    sum_uk = u_hat(:,:,end) + sum_uk - u_hat(:,:,k);
    
    % update first mode's spectrum through wiener filter (on half plane)
    u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
    
    % update first mode's central frequency as spectral center of gravity
    if ~DC
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
    end
    
    % recover full spectrum from analytic signal
    u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
    
    % work on other modes
    for k=2:K
        
        % recompute Hilbert mask
        HilbertMask = (sign(freqs_1*omega(n,1,k) + freqs_2*omega(n,2,k))+1);
        
        % update accumulator
        sum_uk = u_hat(:,:,k-1) + sum_uk - u_hat(:,:,k);
        
        % update signal spectrum
        u_hat(:,:,k) = ((f_hat - sum_uk - mu_hat(:,:)/2).*HilbertMask)./(1+Alpha(k)*((freqs_1 - omega(n,1,k)).^2+(freqs_2 - omega(n,2,k)).^2));
        
        % update signal frequencies
        omega(n+1,1,k) = sum(sum(freqs_1.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        omega(n+1,2,k) = sum(sum(freqs_2.*(abs(u_hat(:,:,k)).^2)))/sum(sum(abs(u_hat(:,:,k)).^2));
        
        % keep omegas on same halfplane
        if omega(n+1,2,k) < 0
            omega(n+1,:,k) = -omega(n+1,:,k);
        end
        
        % recover full spectrum from analytic signal
        u_hat(:,:,k) = fftshift(fft2(real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))))));
    end
    
    % Gradient ascent for augmented Lagrangian
    mu_hat(:,:) = mu_hat(:,:) + tau*(sum(u_hat,3) - f_hat);
    
    % increment iteration counter
    n = n+1;
    
    % convergence?
    uDiff = eps;
    omegaDiff = eps;
    
    for k=1:K
        omegaDiff = omegaDiff + sum(sum(abs(omega(n,:,:) - omega(n-1,:,:)).^2));
        uDiff = uDiff + sum(sum(1/(Hx*Hy)*(u_hat(:,:,k)-u_hat_old(:,:,k)).*conj((u_hat(:,:,k)-u_hat_old(:,:,k)))));
    end
    
    uDiff = abs(uDiff);
    
    u_hat_old = u_hat;

end


%% Signal Reconstruction

% Inverse Fourier Transform to compute (spatial) modes
u = zeros(Hy,Hx,K);
for k=1:K
    u(:,:,k) = real(ifft2(ifftshift(squeeze(u_hat(:,:,k)))));
end;

% Should the omega-history be returned, or just the final results?
%omega = omega(n,:,:);
omega = omega(1:n,:,:);

end

