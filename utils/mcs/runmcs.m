%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% runmcs.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test driver for running MCS 
% with advice on how to choose the tuning parameters
%
% applied to the Jones test functions with default boxes
%
% To solve your own problem, copy this file (to ownmcs.m, say)
% and adapt the first part labelled `problem definition'.
% Then run the driver by typing `ownmcs' at the Matlab prompt.
%
% If you are not satisfied with the results, or if you want to play 
% with the tuning parameters, you also need to adapt the second part
% labelled `change MCS settings'. In particular, for wide bounds,
% you'll probably get better results if you supply your own 
% initialization list.
% 
% On typing `runmcs' at the Matlab prompt,
% the unmodified file produces test results for the six-hump camel
% function; by only changing the data assignment you can also get
% results for the other test functions from Jones et al.
% You may also play with the bounds by modifying the default bounds.
% 

clear; clear mex; clear global; 
format compact;format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define objective function
%
% Each function value f=f(x) is obtained by MCS using a call 
% f=feval(fcn,data,x)
% where data is an arbitrary data object passed to mcs. If you have 
% several data items, put them into a cell array or structure array
% called data.
% If there are no data, use fcn='feval' and the function name as data.
%
fcn = 'feval';
data = 'cam';	% select test function from Jones et al. test set
		% bra  n=2   Branin
		% cam  n=2   six-hump camel
		% gpr  n=2   Goldstein-Price   
		% shu  n=2   Shubert
		% hm3  n=3   Hartman 3
		% s10  n=4   Shekel 10
		% sh5  n=4   Shekel 5
		% sh7  n=4   Shekel 7
		% hm6  n=6   Hartman 6
path(path,'jones');	% add path to directory with function 

% define bounds on variables (+-inf allowed)
%
% u: column vector of lower bounds
% v: column vector of upper bounds
% u(k)<v(k) is required
%
[u,v,fglob] = defaults(data); 	% returns bounds used by Jones et al.
				% and known global optimum
dimension=length(u)		% show dimension
known_global_opt_value=fglob	% show known global minimum value
% u=[-5,0]';v=[10,15]';    		% bra default bounds
% u=[-3,-2]';v=[3,2]';     		% cam default bounds
% u=[-2,-2]';v=[2,2]';     		% gpr default bounds
% u=[-10,-10]';v=[10,10]'; 		% shu default bounds
% u=[0,0,0]';v=[1,1,1]';   		% hm3 default bounds
% u=[0,0,0,0]';v=[10,10,10,10]';   	% sh5,sh7,s10 default bounds
% u=[0,0,0,0,0,0]';v=[1,1,1,1,1,1]';   	% hm6 default bounds
% modify the problem to be unconstrained by activating the next line
% u = -Inf*ones(size(u));v=-u; 


use_defaults=1
% *** If you just want to use the default settings,
% *** you don't need to edit the rest of the file,
% *** If you are not satisfied with the results
% *** (this may happen especially when your box bounds are very wide),
% *** or if you want to play with the tuning parameters,
% *** set use_defaults=0, and modify the rest of the file 
% *** according to your curiosity or ingenuity

if use_defaults, 
  % easy to use version - all parameters preset
  % defaults are being used, it suffices to call  
  %%%%%%%%%%%%%%%%%% simple MCS call %%%%%%%%%%%%%%%%%%
  [xbest,fbest,xmin,fmi,ncall,ncloc]=mcs(fcn,data,u,v);
  % or, with less output,
  % [xbest,fbest]=mcs(fcn,data,u,v);

  xbest	  		% best point found
  fbest     		% best function value
  fglob			% global minimum (known for test functions)
  ncall	  		% number of function values used
  ncloc	  		% number of function values in local searches

  % xmin	  	% columns are points in 'shopping basket'
			% may be good alternative local minima
  % fmi	  		% function values in 'shopping basket'
  nbasket = length(fmi) % number of points in 'shopping basket'

  return;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% change MCS settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flexible use version - all parameters may be modified

% 
% define amount of output printed
prt = 1;	% print level 
		% prt = 0: no output
		% prt = 1: # sweep, minimal nonempty level, # f-calls,
		%          best point and function value
		% prt > 1: in addition levels and function values of
		%          boxes containing the known global minimizers
		%          of a test function

% define global strategy 
%
% smax governs the relative amount of global versus local search. By
% increasing smax, more weight is given to global search.
% 
% Increasing or decreasing stop(1) increases or decreases the amount 
% of work allowed before concluding that nothing more is gained; 
% the default choice is quite conservative, and may try to locate
% a better point long after the global minimum has been found.
% stop(1)=5 works for many easier problems, too, with much fewer
% wasted function values.
% 
% Increasing nf allows to spend more time on boxes that have a chance 
% on their level to contain better points. This may be important for
% hard problems, or for problems with very wide bounds.
% 
% But in each case, it is unclear in advance what change would be most 
% beneficial to a particular problem. 
% We had very mixed experience; if you have many similar problems to 
% solve, the best thing to do is to experiment with a few problems to 
% find the best values, and to use these on the others. 
%
n = length(u);		% problem dimension
smax = 5*n+10;		% number of levels used
nf = 50*n^2; 		% limit on number of f-calls
stop(1) = 3*n;		% = m, integer defining stopping test
stop(2) = -inf;		% = freach, function value to reach
			% if m>0, run until m sweeps without progress
			% if m=0, run until fbest<=freach
			% (or about nf function calls were used)

if 0, 	% known global optimum, for tests only
	% then the entries of stop have a different meaning
  stop(1) = 1.e-4;	% run until this relative error is achieved
  stop(2) = fglob;	% known global optimum value
  stop(3) = 1.e-10;	% stopping tolerance for tiny fglob
end;

% define initialization strategy
%
% for wide boxes, it is advisable (and for unbounded search regions
% strongly advisable) to define a customized initialization list
% that contains for each coordinate at least three reasonable values.
% Without such an initialization list, mcs makes default assumptions
% that roughly amount to estimating reasonable magnitudes from the 
% bounds and in case iinit=1 from assuming order unity if this is 
% within the bounds. 
%
% for a self-defined initialization list, the user should
% write an m-script file init0.m containing a matrix x0 with n
% rows and at least 3 columns and two n-vectors l and L 
% the ith column of x0 contains the initialization list
% values for the ith coordinate, their number is L(i), and
% x0(i,l(i)) is the ith coordinate of the initial point

iinit = 0;	% 0: simple initialization list
		%    (default for finite bounds;
		%     do not use this for very wide bounds)
		% 1: safeguarded initialization list
		%    (default for unbounded search regions)
		% 2: (5*u+v)/6, (u+v)/2, (u+5*v)/6
		% 3: initialization list with line searches
		% else: self-defined initialization list 
		%       stored in file init0.m

% parameters for local search
%
% A tiny gamma (default) gives a quite accurate but in higher 
% dimensions slow local search. Increasing gamma leads to less work 
% per local search but a less accurately localized minimizer
% 
% If it is known that the Hessian is sparse, providing the sparsity 
% pattern saves many function evaluations since the corresponding
% entries in the Hessian need not be estimated. The default pattern
% is a full matrix.
% 
local = 50;		% local = 0: no local search
			% else: maximal number of steps in local search
gamma = eps;		% acceptable relative accuracy for local search
hess = ones(n,n);	% sparsity pattern of Hessian



% defaults are not being used, use the full calling sequence
% (including at least the modified arguments)
%%%%%%%%%%%%%%%%%%%%%%% full MCS call %%%%%%%%%%%%%%%%%%%%%%
[xbest,fbest,xmin,fmi,ncall,ncloc]=...
  mcs(fcn,data,u,v,prt,smax,nf,stop,iinit,local,gamma,hess);

xbest	  		% best point found
fbest     		% best function value
fglob			% global minimum (known for test functions)
ncall	  		% number of function values used
ncloc	  		% number of function values in local searches

% xmin	  		% columns are points in 'shopping basket'
			% may be good alternative local minima
% fmi	  		% function values in 'shopping basket'
nbasket = length(fmi) 	% number of points in 'shopping basket'

