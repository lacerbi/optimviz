function varargout = funlogger(x,varargin)
%FUNLOGGER Record log of calls to a specific function.
%
%   FUNLOGGER([],FUN) initialize logger for function FUN. FUN must be either
%   a string with the function name or a function handle for function F(X). 
%   F takes a numeric array as input and returns a single scalar as output.
%
%   FUNLOGGER([],FUN,N) stores up to N function calls (default N=1e4).
%
%   FUNLOGGER([],FUN,N,OUTSIZE) specifies that FUN returns a numeric array 
%   with OUTSIZE elements (default OUTSIZE=1).
%
%   FVAL = FUNLOGGER(X) returns the value of the logged function at X,
%   where X is a numeric array.
%
%   [FVAL,FUNLOG] = FUNLOGGER(X) also returns the function log structure.
%
%   [...] = FUNLOGGER(X,ARG1,ARG2,...) passes a series of additional 
%   arguments to the call of logged function.
%
%   FUNLOG = FUNLOGGER() returns the function log structure.
%
% 
%   Author: Luigi Acerbi
%   Version: Sep/11/2015
%
persistent funlog;     % Record log of function calls

%% Case 1: No arguments ==> return function log
if nargin < 1
    varargout = {funlog};

    %if isstruct(x)  % Swapped variable order
    %    temp = probstruct;
    %    probstruct = x;
    %    x = temp;
    %end

%% Case 2: Empty X, no output, pass function name or handle ==> Re-initialize log struct
elseif isempty(x) && nargout == 0 && ...
        (ischar(varargin{1}) || isa(varargin{1},'function_handle'))
    funlog = [];
    fun = varargin{1};
    if nargin > 2; N = varargin{2}; else N = []; end
    if nargin > 3; outsize = varargin{3}; else outsize = []; end

    % Store function name and function handle
    if ischar(fun)
        funlog.FuncName = fun;
        funlog.FuncHandle = str2func(fun);
    else
        funlog.FuncName = func2str(fun);
        funlog.FuncHandle = fun;        
    end
    
    % Number of stored function values
    if isempty(N); N = 1e4; end
    
    % Size of function output arguments
    if isempty(outsize); outsize = 1; end
    funlog.OutSize = outsize;
    
    funlog.Clock = tic;     % Time of initialization
    funlog.FuncCount = 0;
    funlog.X = [];
    funlog.Y = [];
    
    funlog.N = N;           % Size of X matrix, temporary
    funlog.last = 0;        % Last entry
    
%% Case 3: Pass X, get output ==> Standard function call    
else
    if isempty(x) && nargin == 2 && isnumeric(varargin{1})
        % MCS call
        x = varargin{1}';
        varargin = [];
        mcs_flag = true;
    else
        mcs_flag = false;
    end
    
    % You need to initialize FUNLOGGER first
    if isempty(funlog)
        error('FUNLOGGER has not been initialized.');
    end

    % Call function
    func = funlog.FuncHandle;

    % if isfield(probstruct,'octstruct'); x = optimct(x,probstruct.octstruct,1); end

    % Check if need to pass probstruct
    try
        if nargin > 1 && ~mcs_flag
            fval = func(x,varargin{:});
        else
            fval = func(x);
        end
    catch funErr
        warning(['Error in executing the logged function ''' funlog.FuncName ''' with input:']);
        x
        rethrow(funErr);
        % fval = NaN(1,funlog.OutSize);
    end

    funlog.FuncCount = funlog.FuncCount + 1;

    % Record log
    if isempty(funlog.X)
        funlog.X = NaN(funlog.N,length(x));
        funlog.Y = NaN(funlog.N,funlog.OutSize);
    end

    funlog.last = max(1,mod(funlog.last+1, funlog.N+1));
    funlog.X(funlog.last,:) = x;
    try
        funlog.Y(funlog.last,:) = fval;
    catch
        pause
    end
    
    varargout{1} = fval;
    if nargout > 1; varargout{2} = funlog; end
end

end