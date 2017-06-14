function [s,stop,alpha,maxmodes,t,method,orderder,liss,sporder] = init(varargin)
%% INIT : internal function for the initialization of the parameters.


s = varargin{1};
N = max(size((s)));
if nargin == 2
  if isstruct(varargin{2})
    inopts = varargin{2};
  else
    error('when using 2 arguments the first one is the analyzed signal and the second one is a struct object describing the options')
  end
elseif nargin > 2
  try
    inopts = struct(varargin{2:end});
  catch
    error('bad argument syntax')
  end
end

% Paramètres par défaut
defopts.stop = 'f';
defopts.alpha = 0.05;
defopts.maxmodes = 8;
defopts.t = 1:max(size(s));
defopts.method = 'os';
defopts.orderder=-1;
defopts.liss = 0;
defopts.sporder = 4;
opt_fields = {'stop','alpha','maxmodes','t','method','orderder','liss','sporder'};
opts = defopts;

if(nargin==1)
  inopts = defopts;
elseif nargin == 0
  error('not enough arguments')
end

% Some checking
names = fieldnames(inopts);
for nom = names'
  if ~any(strcmpi(char(nom), opt_fields))
    error(['bad option field name: ',char(nom)])
  end
  % Et modification des paramètres rentrés
  if ~isempty(eval(['inopts.',char(nom)])) % empty values are discarded
    eval(['opts.',lower(char(nom)),' = inopts.',char(nom),';'])
  end
end

% Mise à jour
stop = opts.stop;
alpha = opts.alpha;
maxmodes = opts.maxmodes;
t = opts.t;
method = opts.method;
orderder = opts.orderder;
liss = opts.liss;
sporder = opts.sporder;

%% Syntax check
% s
if ~isvector(s)
  error('The signal S must have only one row or one column')
end
if size(s,1) > 1
  s = s.';
end

% t
if ~isvector(t)
  error('option field T must have only one row or one column')
end
if ~isreal(t)
  error('time instants T must be a real vector')
end
if size(t,1) > 1
  t = t';
end
if (length(t)~=length(s))
  error('X and option field T must have the same length')
end

% orderder
if max(size(orderder))>1
    if size(orderder,1)>1
        orderder=orderder';
    end
    if size(orderder,2)<maxmodes
        orderder = [orderder orderder(end)*ones(1,maxmodes - length(orderder))];
    end
else
    orderder = orderder(1)*ones(1,maxmodes);
end


% liss
if max(size(liss))>1
    if size(liss,1)>1
        liss=liss';
    end
    if size(liss,2)<maxmodes
        liss = [liss liss(end)*ones(1,maxmodes - length(liss))];
    end
else
    liss = [liss(1) zeros(1,maxmodes - 1)];
end


end
