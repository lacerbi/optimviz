
% Global optimization by multilevel coordinate search (MCS)
% 
% Developed by Waltraud Huyer and Arnold Neumaier
% Dept. of Mathematics, University of Vienna, Austria
%
% Source:
% http://solon.cma.univie.ac.at/~neum/software/mcs/
%
% W. Huyer and A. Neumaier,
% Global optimization by multilevel coordinate search,
% J. Global Optimization 14 (1999), 331-355.
% 
% Please report problems and bugs to
% Arnold Neumaier (neum@cma.univie.ac.at)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCS attempts to find the global minimizer of the bound constrained 
% optimization problem
% 
%  min  f(data,x)
%  s.t. x in [u,v] (a box in R^n),
% 
% where data is a fixed data vector (or other data structure), 
% and f is a function of data and x defined by a user-provided m-File,
% 
% The search is not exhaustive; so the global minimum may be missed.
% However, a comparison to other global optimization algorithms shows
% excellent performance in many cases, especially in low dimensions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before the first use of MCS, read Readme.mcs and runmcs.m; the latter
% is a sample driver for MCS for the test set of Jones et al.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The subdirectory ./jones/ contains the test functions from
% D.R. Jones, C.D. Perttunen and B.E. Stuckman, Journal of Optimization
% Theory and Applications 79 (1993), 157-181
% see jones/Contents.m
%
% the following lists the remaining programs in alphabetical order
% 
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% addloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adds a new point to the list of starting points for local search

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% basket.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a candidate for local search lies in the 'domain of 
% attraction' of a point in the 'shopping basket'
% function [xbest,fbest,xmin,fmi,nbasket,loc,flag] = basket(fcn,data,x,
% f,xmin,fmi,xbest,fbest,stop,nbasket)
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector (or other data structure)
% x(1:n)	candidate for the 'shopping basket' 
% f		its function value
% xmin(1:n,:)  	columns are the base vertices of the boxes in the  
%              	shopping basket
% fmi          	fmi(j) is the function value at xmin(:,j)
% xbest       	current best vertex
% fbest    	current best function value		
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% nbasket	current number of points in the 'shopping basket'
% Output:
% xbest       	current best vertex
% fbest    	current best function value
% xmin(1:n,:)	updated version of the points in the shopping basket
% fmi		their function values
% loc           = 0  candidate lies in the 'domain of attraction' of a
%		     point in the shopping basket 
%		= 1  otherwise
% flag		= 0  the global minimum of a test function has been 
%		     found with the required accuracy relerr
% 		= 1  otherwise
% ncall		number of function calls used in the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% basket1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a candidate for the 'shopping basket' lies in the 
% 'domain of attraction' of a point already in the 'shopping basket'
% function [xbest,fbest,xmin,fmi,loc,flag,ncall] = basket1(fcn,data,x,f,
% xmin,fmi,xbest,fbest,stop,nbasket)
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector (or other data structure)
% x(1:n)	candidate for the shopping basket
% f		its function value
% xmin(1:n,:)  	columns are the base vertices of the boxes in the 
%              	'shopping basket'
% fmi          	fmi(j) is the function value at xmin(:,j)
% xbest       	current best vertex
% fbest    	current best function value
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% nbasket	current number of points in the shopping basket
% Output:
% xbest       	current best vertex
% fbest    	current best function value
% xmin(1:n,:)	updated version of the points in the shopping basket
% fmi		their function values
% loc           = 0  candidate lies in the 'domain of attraction' of a
%		     point in the 'shopping basket 
%		= 1  otherwise
% flag		= 0  the global minimum of a test function has been 
%		     found with the required accuracy prcerr
% 		= 1  otherwise
% ncall		number of function calls used in the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bounds.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [w1,w2] = bounds(n,n0,x,y,u,v)
% computes the bounds for a box with base vertex x and opposite vertex y
% Input:
% n	   dimension of the problem
% n0(1:n)  coordinate i has been split n0(i) times in the history of the
%          box
% x(1:n)   base vertex of the box
% y(1:n)   opposite vertex of the box
% [u,v]    original box
% Output:
% w1(1:n)  vector of lower bounds
% w2(1:n)  vector of upper bounds
% i.e. the box is [w1,w2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chkloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a point has already been used as starting point for a
% local search; in that case loc is set to 0, otherwise, it is set to 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chrelerr.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function flag = chrelerr(fbest,stop)
% checks whether the required tolerance for a test function with known
% global minimum has already been achieved
% Input:
% fbest		function value to be checked
% stop(1)	relative error with which a global minimum with not too
%		small absolute value should be reached
% stop(2)	global minimum function value of a test function
% stop(3)	if abs(fglob) is very small, we stop if the function
%		value is less than stop(3)
% Output:
% flag          = 0 the required tolerance has been achieved
% 		= 1 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chvtr.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function chvtr(f,vtr)
% checks whether a required value to reach has already been reached; in
% that case flag is set to 0, otherwise it is set to 1
% Input:
% f	function value to be checked
% vtr	value to reach
% Output:
% flag	= 0  vtr has been reached
%	= 1  otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% csearch.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xmin,fmi,g,G,nfcsearch] = csearch(fcn,data,x,f,u,v,hess)
% makes local line searches with a maximally 6 points in all coordinate
% directions
% Input:
% fcn = 'fun'	name of function fun(data,x), x an n-vector
% data		data vector
% x		starting point for the coordinate searches
% f		its function value
% [u,v]		box bounds for the optimization problem
% hess		sparsity pattern of the Hessian (default: hess = 
%		ones(n,n), 
%		where n is the dimension of the problem)
%
% Output:
% xmin		end point after the coordinate searches
% fmi		its function value
% g		estimated gradient at xmin
% G		estimated Hessian at xmin
% nfcsearch	number of function values used in the coordinate 
%		searches
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% exgain.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [e,isplit,splval] = exgain(n,n0,l,L,x,y,x1,x2,fx,f0,f1,f2)
% determines the splitting index, the splitting value and the expected
% gain vector e for (potentially) splitting a box by expected gain
% Input:
% n        dimension of the problem
% n0(1:n)  the ith coordinate has been split n0(i) times in the history
%          of the box
% l(1:n)   pointer to the initial point of the initialization list
% L(1:n)   lengths of the initialization list
% x(1:n)   base vertex of the box
% y(1:n)   opposite vertex of the box
% x1(1:n), x2(1:n), f1(1:n), f2(1:n)
%          x1(i) and x2(i) and the corresponding function values f1(i)
%          and f2(i) used for quadratic interpolation in the ith
%          coordinate 
% fx       function value at the base vertex
% f0(1:max(L),1:n)  function values appertaining to the init. list
% Output:
% e(1:n)   e(i) maximal expected gain in function value by changing 
%          coordinate i
% isplit   splitting index
% splval   = Inf  if n0(isplit) = 0
%          = splitting value  otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fbestloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a new point in the shopping basket (found by local 
% search) is better than the current best point; in that case it 
% updates the best point and its function value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% genbox.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ipar,level,ichild,f] = genbox(par,level0,nchild,f0)
% generates a box with parent box # ipar = par, level = level0, ichild 
% = nchild and base vertex function value f = f0 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hessian.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h = hessian(i,k,x,x0,f,f0,g,G)
% computes the element G(i,k) of the Hessian of the local quadratic 
% model
% Input:
% i, k		indices (k < i)
% x(1:n)	'neighbor' of x0 used for computing G(i,k); differs from
%		x0 only in the ith and kth component
% x0(1:n)	point around which the quadratic model is computed
% f		function value at x
% f0		function value at x0
% g(1:i)	components of the gradient of the local quadratic model %		that have already been computed
% G		components of the Hessian of the local quadratic model
%		that have already been computed
% Output:
% h = G(i,k)	newly computed nondiagonal element of the Hessian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% init.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f0,istar,ncall] = init(fcn,data,x0,l,L,n)
% computes the function values corresponding to the initialization list
% and the pointer istar to the final best point x^* of the init. list
% Input:
% fcn = 'fun' 	name function fun(data,x), x an n-vector
% data	      	data vector 
% x0(1:n,1:max(L))  coordinates used in the initialization list  
% l(1:n)      	pointer to the initial point x^0 = (x0(1,l(1)),...,
%             	x0(n,l(n)) of the initialization list
% L(1:n)      	the init. list contains L(i) values for coordinate i
% n           	dimension of the problem
% Output:
% f0(1:max(L),1:n)  function values appertaining to the initialization
%             	list; f0(j,i) = function value at x = (x0(1,istar(1)),
%		...,x0(i-1,istar(i-1)),x0(i,j),x0(i+1,l(i+1)),...,
%		x0(n,l(n)))
% istar(1:n)  	pointer to the final best point x^* of the init. list
%             	x^* = (x0(1,istar(1)),...,x0(n,istar(n)))
% ncall	      	number of function calls used in the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% init0.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% self-defined initialization list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initbox.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ipar,level,ichild,f,isplit,p,xbest,fbest]=initbox(x0,f0,l,
% L,istar,u,v,prt)
% generates the boxes in the initialization procedure
% For each coordinate i we compute the quadratic interpolant to any
% three consecutive (x0(i,j),f0(j,i)). xl and xu are the minimum and
% maximum of the quadratic interpolant, respectively. If x0(i,istar(i)) 
% belongs to two boxes, the one containing the minimizer of the 
% quadratic model in coordinate i is taken as the current box for the
% next coordinate. All boxes generated in the initialization procedure
% get the level of their parent box increased by one.
%
% Input:
% x0(1:n,1:max(L))   coordinates used in the initialization list
% f0(1:max(L),1:n)  function values appertaining to the initialization
%             list; f0(j,i) = function value at x = (x0(1,istar(1)),...,
%             x0(i-1,istar(i-1)),x0(i,j),x0(i+1,l(i+1)),...,x0(n,l(n)))
% l(1:n)      pointer to the initial point x^0 = (x0(1,l(1)),...,
%             x0(n,l(n)) of the initialization list
% L(1:n)      the init. list contains L(i) values for coordinate i
% istar(1:n)  pointer to the final best point x^* of the init. list
%             x^* = (x0(1,istar(1)),...,x0(n,istar(n))) 
% [u,v]       initial box
% Output:
% ipar        vector containing the labels of the parents of the newly
%             created boxes
% level       vector containing their levels
% ichild      the absolute values of this vector specify which child 
%             they are; ichild(j) < 0 for the boxes generated in the 
%             initialization procedure
% f(1:2,:)    f(1,j) contains the base vertex function value of box j
% isplit      vector containing the splitting indices of the boxes that
%             have already been split
% p(1:n)      ranking of variability; the ith component has the p(i)th 
%             highest estimated variability 
% xbest       best vertex after the initialization procedure
% fbest       best function value after the initialization procedure
% prt	      print level

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initlist.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x0,f0,l,L,istar,ncall] = initlist(fcn,data,u,v)
% generates an initialization list with the aid of line searches
%
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector
% [u,v]       	box in which the optimization is carried out (u, v 
%             	n-vectors)
% Output:
% x0		array with n rows and at least 3 columns; the ith row
%		contains the initialization list values for the ith 
%		coordinate 
% f0		corresponding function values
% l		x0(i,l(i)) is the ith coordinate of the initial point
% L		
% istar		x0(i,istar(i)) is the ith coordinate of the final best 
%		point
% ncall		number of function calls used in the program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lsearch.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xmin,fmi,ncall,flag] = lsearch(fcn,data,x,f,f0,u,v,nf,stop,
% maxstep,gamma,hess)
% Input:
% fcn = 'fun'	name of function fun(data,x), x an n-vector
% data		data vector
% x		starting point for the local search
% f		its function value
% f0		an additional 'typical' function value (f0 >= f)
% [u,v]		original box
% nf		limit on function calls
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% maxstep	maximum number of steps in local search (default: 50)
% gamma		stopping criterion for local search (default: eps)
%		the local search is stopped if abs(g)'*max(abs(x),
%		abs(xold)) < gamma*(f0-f)
% hess		sparsity pattern of the Hessian (default hess = 
%		ones(n,n), where n is the dimension of the problem)
%
% Output:
% xmin(1:n)	result of the local search
% fmi		its function value
% ncall		number of function evaluations
% flag      	is set to 0 if the global minimum of a test function 
%		is found with the desired accuracy; otherwise it is 
%		set to one

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mcs.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xbest,fbest,xmin,fmi,ncall,ncloc,flag]=mcs(fcn,data,u,v,prt,
% smax,nf,stop,iinit,local,gamma,hess)
% MCS global optimization for function defined by fcn in the 
% n-dimensional box [u,v]
%
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector (or other data structure)
% [u,v]       	box in which the optimization is carried out (u, v 
%             	n-vectors)
% prt		print level
% 		prt = 0: no printing
% 		prt = 1: # sweep, minimal nonempty level, # f-calls, 
% 		best point and function value (default)
% 		prt > 1: only meaningful for test functions with known
% 		global minimizers
% 		in addition levels and function values of boxes 
% 		containing the global minimizers of a test function
% smax        	number of levels (default: 5*n+10)
% nf         	maximum number of function evaluations (default: 50*n^2)
% stop         	stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
%            	(default: stop = 3*n)
% iinit       	parameter defining the initialization list
%             	= 0        corners and midpoint (default for finite u,v)
%             	= 1        safeguarded version *default otherwise)
% 		= 2        5u/6 + v/6, u/6 + 5v/6 and midpoint
%		= 3        initialization list with line searches
%             	otherwise  self-defined init. list 
%		for a self-defined initialization list, the user should
% 		provide an m-script file init0.m containing a matrix x0 
%		with n rows and at least 3 columns and two n-vectors l 
%		and L 
%		the ith column of x0 contains the initialization list
%		values for the ith coordinate, their number is L(i), and
%		x0(i,l(i)) is the ith coordinate of the initial point
% local		local = 0: no local search
%		otherwise: maximal number of steps in local search
%		(default: 50) 
% gamma		stopping criterion for local search (default: eps)
%           	the local search is stopped if abs(g)'*max(abs(x),
%		abs(xold)) < gamma*(f0-f) 
% hess		sparsity pattern of the Hessian for local search 
%		(default: hessian = ones(n,n))
%
% Output:
% xbest(1:n)  	current best point 
% fbest    	function value at xbest
% xmin        	matrix with n rows; the columns are the points in the
%             	'shopping basket' (i.e. good points resp. local 
%		minimizers)
% fmi         	function values corresponding to the 'shopping basket';
%		fmi(i) is the function value at xmin(:,i)
% ncall       	number of function evaluations
% ncloc		number of function evaluations used for local search
% flag        	specifies which stopping criterion has been used
%             	= 0  a (known) global minimum fglob of a test function 
%                    has been found with the required relative error 
%		     relerr
%             	= 1  the division procedure has been completed
%             	= 2  the maximum number nf of function calls has been
%                    reached without finding a known minimum with the
%                    required relative error or completing the division
%                    procedure
%		= 3  stop(1) sweeps without progress (for stop(1) >= 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% neighbor.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2] = neighbor(x,delta,u,v)
% computes 'neighbors' x1 and x2 of x needed for making triple search 
% and building a local quadratic model such that x(i), x1(i), x2(i) are
% pairwise distinct for i = 1,...,n
% Input:
% x  	points for which 'neighbors' should be found
% delta	distance of the neighbors
% [u,v]	original box for global optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% polint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function d = polint(x,f)
% quadratic polynomial interpolation
% Input: 
% x(1:3)  3 pairwise distinct support points
% f(1:3)  corresponding function values
% Output:
% d(1:3)  the interpolating polynomial is given by
%         p(x) = d(1) + d(2)(x - x(1)) + d(3)(x - x(1))(x - x(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% polint1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [g,G] = polint1(x,f)
% quadratic polynomial interpolation
% Input:
% x(1:3)  3 pairwise distinct support points
% f(1:3)  corresponding function values
% Output:
% g, G
% the interpolating polynomial is given by
% p(x) = f(1) + g(x - x(1)) + G/2(x - x(1))^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quadmin.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x = quadmin(a,b,d,x0)
% computes the minimum x of the quadratic polynomial
% p(x) = d(1) + d(2)(x - x0(1)) + d(3)(x - x0(1))(x - x0(2))
% in the interval [a,b]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quadpol.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f = quadpol(x,d,x0)
% evaluates the quadratic polynomial
% p(x) = d(1) + d(2)(x - x0(1)) + d(3)(x - x0(1))(x - x0(2))
% at x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% range.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function range(u,v,x,p)
% computes the range [amin,amax] for making a line search along x + ap
% when x + ap is restricted to [u,v]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% runmcs.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test file for running MCS 
% with the Jones test functions with default boxes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splinit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xbest,fbest,f0,xmin,fmi,ipar,level,ichild,f,flag,ncall] = 
% splinit(fcn,data,i,s,smax,par,x0,n0,u,v,x,y,x1,x2,L,l,xmin,fmi,xbest,
% fbest,ipar,level,ichild,f,stop,prt)
% splits box # par at level s according to the initialization list
% in the ith coordinate and inserts its children and their parameters
% in the list 
% Input:
% fcn = 'fun'  	name of function fun(data,x), x an n-vector
% data		data vector
% i            	splitting index
% s            	level of the box
% smax         	depth of search
% par          	label of the box
% x0(1:n,1:max(L)) initialization list
% n0(1:n)      	coordinate i has been split n0(i) times in the history
%              	of the box
% [u,v]        	original box
% x(1:n)       	base vertex of the box
% y(1:n)       	opposite vertex
% x1(1:n), x2(1:n) 'neighbors' of x such that x(i), x1(i), x2(i) are
%		pairwise distinct for i = 1,...,n
% L(1:n)       	lengths of the initialization list
% l(1:n)       	pointer to the initial point in the initialization list
% xmin(1:n,:)  	columns are the base vertices of the boxes in the  
%              	'shopping basket'
% fmi          	fmi(j) is the function value at xmin(:,j)
% xbest       	current best vertex
% fbest    	current best function value
% ipar         	vector containing the labels of the parents of the boxes
%              	not in the shopping basket
% level        	vector containing their levels
% ichild       	the absolute values of this vector specify which child 
%              	they are; ichild(j) < 0 if box j was generated by 
%              	splitting according to the init. list (in the init.
%              	procedure or later) and ichild(j) > 0 otherwise
% f(1:2,:)     	f(1,j) is the base vertex function value of box j and
%              	f(2,j) contains the function value at its splitting 
%              	value (if box j has already been split by default)	
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% prt		print level
% Output:
% xbest       	current best vertex
% fbest    	current best function value
% f0(1:L(i))   	base vertex function values of the newly created boxes 
% xmin(1:n,:)  	as before plus newly created boxes 
% fmi          	as before plus newly created boxes 
% ipar         	as before plus newly created boxes
% level        	as before plus newly created boxes
% ichild       	as before plus newly created boxes
% f            	as before plus newly created boxes
% flag         	output flag
%              	= 0 if the known global minimum of a test function has 
%                   been found with the required relative error
%              	= 1 otherwise 
% ncall		number of function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xbest,fbest,xmin,fmi,ipar,level,ichild,f,flag,ncall] = 
% split(fcn,data,i,s,smax,par,n0,u,v,x,y,x1,x2,z,xmin,fmi,ipar,level,
% ichild,f,xbest,fbest,stop,prt)
% splits box # par at level s in its ith coordinate into two or three
% children and inserts its children and their parameters in the list
%
% Input:
% fcn = 'fun'  	name of function fun(data,x), x an n-vector
% data		data vector
% i            	splitting index
% s            	level of the box to be split
% smax         	number of levels
% par          	label of the box
% n0(1:n)      	coordinate i has been split n0(i) times in the history
%              	of the box
% [u,v]        	original box
% x(1:n)       	base vertex of the box
% y(1:n)       	opposite vertex
% x1(1:n), x2(1:n) 'neighbors' of x such that x(i), x1(i), x2(i) are
%		pairwise distinct for i = 1,...,n
% z(1:2)       	z(1) = value of the ith coordinate of the base vertex
%              	z(2) = splitting value
%              	z(2) ~= y(i) split into 3 children
%              	z(2) = y(i)  split into 2 children
% xmin(1:n,:)  	columns are the base vertices of the boxes in the base 
%              	list
% fmi          	fmi(j) is the function value at xmin(:,j)
% ipar         	vector containing the labels of the parents of the boxes
%              	not in the `shopping basket'
% level        	vector containing their levels
% ichild       	the absolute values of this vector specify which child 
%              	they are; ichild(j) < 0 if box j was generated by 
%              	splitting according to the init. list (in the init.
%              	procedure or later) and ichild(j) > 0 otherwise
% f            	f(1,j) is the base vertex function value of box j and
%              	f(2,j) contains the function value at its splitting 
%              	value (if box j has already been split by default)
% xbest(1:n)  	current best vertex
% fbest    	current best function value
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% prt		print level
% Output:
% xbest(1:n) 	current best vertex
% fbest    	current best function value
% xmin(1:n,:)  	as before plus newly created boxes
% fmi          	as before plus newly created boxes 
% ipar         	as before plus newly created boxes
% level        	as before plus newly created boxes
% ichild       	as before plus newly created boxes
% f            	as before plus newly created boxes; in particular 
%              	f(2,par) and f(1,:) for the newly created boxes are
%              	defined in this program
% flag         	output flag
%              	= 0 if the known global minimum of a test function has 
%                   been found with the required relative error
%               = 1 otherwise  
% ncall		number of function evaluations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function z = split1(x1,x2,f1,f2)
% Input: points x1 and x2, x1 < x2, and corresponding function values f1
%        and f2
% splits the interval [x1,x2] according to the golden section rule
% the part containing the better point gets the larger fraction of the 
% interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x1 = split2(x,y)
% determines a value x1 for splitting the interval [min(x,y),max(x,y)]
% is modeled on the function subint with safeguards for infinite y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% splrnk.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [isplit,splval] = splrnk(n,n0,p,x,y)
% determines the splitting index and splitting value for splitting a
% box by rank
% Input:
% n        dimension of the problem
% n0(1:n)  coordinate i has been split n0(i) times in the history of the
%          box to split
% p(1:n)   ranking of estimated variability of the function in the 
%          different coordinates
% x(1:n)   base vertex of the box
% y(1:n)   opposite vertex of the box
% Output:
% isplit   splitting index
% splval   = Inf  if n0(isplit) = 0 (indicates that the box has to be
%                 split according to the initialization list)
%          = splitting value  otherwise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strtsw.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function s = strtsw(smax,level,f)
% updates the record list for starting a new sweep and computes the
% lowest level containing non-split boxes

% Input:
% smax         depth of search
% level(1:nboxes)  levels of these boxes (level(j) = 0 if box j has 
%              already been split)
% f(1:nboxes)  their function values at the base vertices

% Output:
% s            lowest level with record(s) ~= 0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2] = subint(x,y)
% computes for real x and real or infinite y two points x1 and x2 in
% [min(x,y),max(x,y)] that are neither too close nor too far away from x


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% triple.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xtrip,ftrip,g,G,x1,x2,nf] = triple(fcn,data,x,f,x1,x2,u,v,G)
% finding a local quadratic model by triple search
% Input:
% fcn = 'fun'	name of function fun(data,x), x an n-vector
% x(1:n)	start vector
% f		its function value
% x1(1:n), x2(1:n)  'neighbors' of x with x(i), x1(i), x2(i) pairwise
%		distinct for i = 1,...,n
% [u,v]		original box
% G		Hessian of the quadratic model around xtrip
% Output:
% xtrip(1:n)	best point found 
% ftrip		its function value
% g, G		gradient and Hessian of the quadratic model around xtrip
% the quadratic model is given by
% q(x) = ftrip + g^T(x-xtrip) + 0.5(x-xtrip)G(x-xtrip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f1,f2,fold] = updtf(n,i,x1,x2,f1,f2,fold,f)
% subprogram used by vertex.m
% If we leave the i-coordinate x(i) of the base vertex of a box by
% going back to the base vertices of its ancestors, we have to 
% generate 'fictitious' function values by assuming a separable model
% for all coordinates i1 ~= i for which coordinate and function values
% have not yet been found
% Input:
% n        dimension of the problem
% i        splitting index
% x1(1:n), x2(1:n), f1(1:n), f2(1:n)
%          x1(j) and x2(j) and the corresponding function values f1(j)
%          and f2(j) used for quadratic interpolation in the jth
%          coordinate 
%          x1(j) = Inf or x2(j) = Inf indicates that these quantities
%          have not been found yet
% fold     base vertex function value of the previously considered box
% f        base vertex function value of the current box
% Output
% f1(1:n), f2(1:n) updates of f1 and f2 
% fold     base vertex function value of the current box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtoptl.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function updtoptl(i,x,y,iopt,level,f)
% a daughter box of a box containing (a) local minimizer(s) is checked
% whether it contains any minimizers, and in case it contains any, the
% corresponding levels and function values are updated 
% Input:
% i      splitting index of the parent box
% x,y    bounds of the ith coordinate of the daughter box 
% iopt   vector containing the indices of the global minimizers the
%        parent box contained
% level  level of the daughter box
% f      base function value of the daughter box

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtrec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function updtrec(j,s,f)
% Input:
% j           label of a box
% s           its level
% f           vector containing the base vertex function values of the
%             already defined boxes
% updates the pointer record(s) to the best non-split box at level s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vert1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,x1,x2,f1,f2] = vert1(j,z,f,x1,x2,f1,f2)
% subprogram used by vertex.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vert2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2,f1,f2] = vert2(j,x,z,f,x1,x2,f1,f2)
% subprogram used by vertex.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vert3.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2,f1,f2] = vert3(j,x0,f0,L,x1,x2,f1,f2)
% subprogram called by vertex.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vertex.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [n0,x,y,x1,x2,f1,f2] = vertex(j,n,u,v,v1,x0,f0,ipar,isplit,
% ichild,z,f,l,L)
% computes the base vertex x and the opposite vertex y of the box # j 
% of MCS and the 'neighboring vertices' x1 and x2 and their function
% values f1 and f2 needed for separable quadratic interpolation
% Input:
% j       label of the box whose vertices etc. are to be computed 
%         (in the main program mcs.m)
% n       dimension of the problem
% [u,v]   original box
% v1      opposite vertex of the original box
% x0(1:n,1:max(L)) coordinates used in the initialization list
%         x^0 = (x0(1,l(1)),...,x0(n,l(n)) base vertex of the original 
%	  box
% f0(1:max(L),1:n) function values appertaining to the initialization 
%	  list
%         f0(m,i) is the function value at x = (x0(1,istar(1)),..., 
%         x0(i-1,istar(i-1)),x0(i,m),x0(i+1,l(i+1)),...,x0(n,l(n)))
% f0(1:max(L),m), m > n  function values appertaining to a later box
%         split according to the initialization list
% ipar(m) label of the parent of box m
% isplit(m) = -(splitting index) of box m  if it is split in the
%             initialization procedure
%           = splitting index of box m  otherwise
% ichild(m) = -(number which child box m is)  if box m was generated
%             by splitting according to the initialization list (in
%             the initialization procedure or later)
%           = number which child box m is  otherwise
% z(1:2,:)  z(1,m) = value of the isplit(m)th coordinate of the base
%           vertex of box m
%           z(2,m) = k  if box m is split according to the 
%                       initialization list and f0(:,k) contains the
%                       function values obtained by splitting the box
%           z(2,m) = splitting value of box m  otherwise
% f(1:2,:)  f(1,m) = function value at the base vertex of box m
%           f(2,m) = function value at the splitting point of box m
%                    (if z(2,m) ~= Inf)
% l(1:n)  pointer defining the initial base vertex x^0 (see above)
% L(1:n)  the init. list contains L(i) values for coordinate i
% Output:
% n0(1:n)  n0(i) indicates that the ith coordinate has been split n0(i)
%          times in the history of box j
% x(1:n)   base vertex of box j
% y(1:n)   opposite vertex of box j
% x1(1:n), x2(1:n), f1(1:n), f2(1:n)
%          x1(i) and x2(i) and the corresponding function values f1(i)
%          and f2(i) are used for quadratic interpolation in the ith
%          coordinate 
% Since it is not always possible to find for each i two vertices that 
% differ from x only in coordinate i, we have to generate 'fictitious'
% vertices and function values using the assumption that the function
% is separable


