function [Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin)
% Function brevridge_mult : extracts the ridges of a multicomponent signal
% Inputs:
%   Tx SQ transform of s
%   lfs : log(fs) frequencies
%   nr : number of ridges
%   lambda : lambda parameter
%   clwin : frequency clearing window
% Outputs:
%   Cs Cs(:,j) is the j-th ridge location

if nargin<5
    lambda=0.001;
    clwin=1;
    doplot=0;
elseif nargin<6
    clwin=1;
    doplot=0;
elseif nargin<7
    doplot=0;
end

Txs = Tx;

[na,N] = size(Tx);

Cs = zeros(N, nr);
Es = zeros(nr, 1);

for j=1:nr
    [Cs(:,j), Es(j)] = brevridge(Tx, fs, lambda);

    % Remove this curve from the representation
    % Max/min frequencies for each time step
    for cli=0:clwin
        fb = max(1, Cs(:,j)-cli);
        fe = min(na, Cs(:,j)+cli);
        Tx([0:N-1]'*na+fb) = sqrt(eps);
        if cli>0, Tx([0:N-1]'*na+fe) = sqrt(eps); end
    end
end

Cs = Cs';

end
