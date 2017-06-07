function [x,fval,outflag] = cmaes_wrapper(fun,x0,LB,UB,PLB,PUB,options,varargin)
%CMAES_WRAPPER A wrapper for the CMA-ES algorithm
%
%  X = CMAES_WRAPPER(FUN,X0,LB,UB) returns the minimum X of a function with 
%  function handle FUN. LB and UB are lower/upper bounds.
%
%  X = CMAES_WRAPPER(FUN,X0,LB,UB,PLB,PUB) uses PLB and PUB as plausible
%  lower and upper bounds (used for initialization). If PLB (resp. PUB) are
%  empty, uses LB (resp. UB) instead.
%
%  X = CMAES_WRAPPER(FUN,X0,LB,UB,PLB,PUB,OPTIONS) passes structure
%  OPTIONS (see 'cmaes' for details).
%
%  X = CMAES_WRAPPER(...,DATA) additional arguments are passed to FUN.
%
%  [X,FVAL] = CMAES_WRAPPER(...) also returns the value of the function at
%  the minimum.
%
%  See also CMAES.

%  Author: Luigi Acerbi
%  Email: luigi.acerbi@gmail.com
%  Date: May 18, 2016

if nargin < 5; PLB = []; end
if nargin < 6; PUB = []; end
if nargin < 7; options = []; end

if isempty(PLB); PLB = LB; end
if isempty(PUB); PUB = UB; end

nvars = numel(x0);

% Set maximum number of function evaluations
if ~isfield(options,'MaxFunEvals') || isempty(options.MaxFunEvals)
    options.MaxFunEvals = 1000*nvars;
end

options.Lbounds      = LB(:);  % Lower bound
options.Ubounds      = UB(:);  % Upper bound
options.Restarts     = 1;

% Disable all writing
options.ReadSignals  = 'off';
options.SaveVariables = 'off';
options.LogModulo    = 0;
options.LogTime      = 0;

% Manage verbosity
if isfield(options,'Display') && (strcmpi(options.Display,'iter') || strcmpi(options.Display,'on'))
    options.DispFinal  = 'on';
    options.DispModulo = ceil(options.MaxFunEvals/100);
else
    options.DispFinal  = 'off';
    options.DispModulo = Inf;
end
if isfield(options,'Display'); options = rmfield(options,'Display'); end

% Use active covariance reduction (almost always more efficient)
options.CMA.Active = 1;

% Suggested starting sigma (std of uniform distribution over the range)
if ~isempty(PLB) && ~isempty(PUB)
    insigma = (PUB(:) - PLB(:))/sqrt(12);
else
    insigma = 10;
end

% if ~ischar(fun); strfun = func2str(fun); else strfun = fun; end

if ischar(fun); fun = str2func(fun); end
cmaes_fun(fun);

[~,~,counteval,outflag,out] = ...
    cmaes('cmaes_fun',x0',insigma,options,varargin{:});

% Take best solution ever found
x = out.solutions.bestever.x;
fval = out.solutions.bestever.f;

if size(x0,1) > 1; x = x'; end