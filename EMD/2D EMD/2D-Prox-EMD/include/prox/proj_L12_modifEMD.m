function wp2 = proj_L12_modifEMD(y,eta)

%function y = project_tv(x, eta, w)
%
%  Created on: --/--/--, Mireille El Gheche
% Modified on: 19/01/12, Giovanni Chierchia
% Modified on: 10/02/12, Nelly Pustelnik
%
% The function computes the projection over the TV constraint.

[n,m,~] = size(y);
x = zeros(n*m,2);
x(:,1) = reshape(y(:,:,1),n*m,1);
x(:,2) = reshape(y(:,:,2),n*m,1);

w=[];


% default input
if nargin < 3
    w = [];
end


% unpack the gradients
[N, Nblock] = size(x);
x_c = x;

% compute the actual TV
mod = sum( sqrt(sum(x_c.^2, 2)) );

% make the projection (only if mod > eta)
if mod <= eta
    
    y = x_c;
    
else
    % compute the projection
    xa  = sqrt( sum(x_c.^2, 2) );
    idx = xa < eps;
    ya  = oneProjector(xa, w(:), eta);
    ya  = ya ./ xa; 
    %ya(idx) = 0; %temporaire
    y = x_c .* repmat(ya, [1 Nblock]);
end
wp = y;

wp2 = zeros(n,m,2);
wp2(:,:,1) = reshape(wp(:,1),n,m);
wp2(:,:,2) = reshape(wp(:,2),n,m);


% ----------------------------------------------------------------------
function [x,itn] = oneProjector(b,d,tau)
% ----------------------------------------------------------------------
% ONEPROJECTOR  Projects b onto the weighted one-norm ball of radius tau
%
%    [X,ITN] = ONEPROJECTOR(B,TAU) returns the orthogonal projection
%    of the vector b onto the one-norm ball of radius tau. The return
%    vector X which solves the problem
%
%            minimize  ||b-x||_2  st  ||x||_1 <= tau.
%               x
%
%    [X,ITN] = ONEPROJECTOR(B,D,TAU) returns the orthogonal
%    projection of the vector b onto the weighted one-norm ball of
%    radius tau, which solves the problem
%
%            minimize  ||b-x||_2  st  || Dx ||_1 <= tau.
%               x
%
%    If D is empty, all weights are set to one, i.e., D = I.
%
%    In both cases, the return value ITN given the number of elements
%    of B that were thresholded.
%
% See also spgl1.


% Check arguments
if nargin < 2
  error('The oneProjector function requires at least two parameters');
end
if nargin < 3
  tau = d;
  d   = [];
end

% Check weight vector
if isempty(d), d = 1; end;

if ~isscalar(d) && ( length(b) ~= length(d) )
  error('Vectors b and d must have the same length');
end

% Quick return for the easy cases.
if isscalar(d)  &&  d == 0
   x   = b;
   itn = 0;
   return
end

% Get sign of b and set to absolute values
s = sign(b);
b = abs(b);

% Perform the projection
if isscalar(d)
  [x,itn] = oneProjectorMex(b,tau/d);
else
  d   = abs(d);
  idx = find(d > eps); % Get index of all non-zero entries of d
  x   = b;             % Ensure x_i = b_i for all i not in index set idx
  [x(idx),itn] = oneProjectorMex(b(idx),d(idx),tau);
end

% Restore signs in x
x = x.*s;




% ----------------------------------------------------------------------
function [x, itn] = oneProjectorMex(b,d,tau)
% ----------------------------------------------------------------------
%function [x, itn] = oneProjectorMex(b,d,tau) 
% Return the orthogonal projection of the vector b >=0 onto the
% (weighted) L1 ball. In case vector d is specified, matrix D is
% defined as diag(d), otherwise the identity matrix is used.
%
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||Dx||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1, oneProjector.

if nargin < 3
    tau = d;
    d   = 1;
end

if isscalar(d)
    [x,itn] = oneProjectorMex_I(b,tau/abs(d));
else
    [x,itn] = oneProjectorMex_D(b,d,tau);
end




% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I(b,tau)
% ----------------------------------------------------------------------

% Initialization
n     = length(b);
x     = zeros(n,1);
bNorm = norm(b,1);

% Check for quick exit.
if (tau >= bNorm), x = b; itn = 0; return; end
if (tau <  eps  ),        itn = 0; return; end

% Preprocessing (b is assumed to be >= 0)
[b,idx] = sort(b,'descend'); % Descending.

csb       = -tau;
alphaPrev = 0;
for j= 1:n
    csb       = csb + b(j);
    alpha     = csb / j;
    
    % We are done as soon as the constraint can be satisfied
    % without exceeding the current minimum value of b
    if alpha >= b(j)
        break;
    end
    
    alphaPrev = alpha;
end

% Set the solution by applying soft-thresholding with
% the previous value of alpha
x(idx) = max(0,b - alphaPrev);

% Set number of iterations
itn = j;




% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_D(b,d,tau)
% ----------------------------------------------------------------------

% Initialization
n = length(b);
x = zeros(n,1);

% Check for quick exit.
if (tau >= norm(d.*b,1)), x = b; itn = 0; return; end
if (tau <  eps         ),        itn = 0; return; end

% Preprocessing (b is assumed to be >= 0)
[bd,idx] = sort(b ./ d,'descend'); % Descending.
b  = b(idx);
d  = d(idx);

% Optimize
csdb = 0; csd2 = 0;
soft = 0; alpha1 = 0; i = 1;
while (i <= n)
    csdb = csdb + d(i).*b(i);
    csd2 = csd2 + d(i).*d(i);
    
    alpha1 = (csdb - tau) / csd2;
    alpha2 = bd(i);
    
    if alpha1 >= alpha2
        break;
    end
    
    soft = alpha1;  i = i + 1;
end
x(idx(1:i-1)) = b(1:i-1) - d(1:i-1) * max(0,soft);

% Set number of iterations
itn = i;






