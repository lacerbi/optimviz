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
%             	(default: stop = 3*n)
% iinit       	parameter defining the initialization list
%             	= 0        corners and midpoint (default for finite u,v)
%             	= 1        safeguarded version *default otherwise)
% 		= 2        5u/6 + v/6, u/6 + 5v/6 and midpoint
%		= 3        initialization list with line searches
%             	otherwise  self-defined init. list (to be stored in 
%			   init0.m)
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
%
% Uses the following m-files (directly or indirectly):
% addloc.m 
% basket.m
% basket1.m
% bounds.m
% chkloc.m
% chrelerr.m
% chvtr.m
% csearch.m 	called by lsearch.m
% exgain.m
% fbestloc.m
% genbox.m
% hessian.m	called by triple.m
% init.m
% initbox.m
% initlist.m
% gls.m and its subprograms   called by lsearch.m 
% lsearch.m
% minq.m and its subprograms  called by lsearch.m
% neighbor.m
% polint.m	called by exgain.m and initbox.m
% polint1.m	called by triple.m
% quadmin.m
% quadpol.m
% range.m	called by lsearch.m
% splinit.m
% split.m
% split1.m
% split2.m   	called by splrnk.m
% splrnk.m
% strtsw.m
% subint.m
% triple.m	called by lsearch.m
% updtf.m    	called by vertex.m
% updtoptl.m
% updtrec.m   	called by splinit.m and split.m
% vert1.m    	called by vertex.m
% vert2.m    	called by vertex.m
% vert3.m    	called by vertex.m
% vertex.m




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xbest,fbest,xmin,fmi,ncall,ncloc,flag]=mymcs(fcn,data,u,v,prt,smax,nf,stop,iinit,local,gamma,hess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global variables
global foptbox nbasket nboxes ncall nglob nsweep nsweepbest optlevel  record xglob xloc
% foptbox(1:nglob)  function value(s) of the box(es) containing the (a)
%             	global minimizer of a test function
% nbasket   	counter for boxes in the 'shopping basket'
% nboxes      	counter for boxes not in the 'shopping basket'
% nglob       	number of global minimizers of a test function
% nloc		(for local ~= 0) counter of points that have been used
% 		as starting points for a local search
% nsweep      	sweep counter
% nsweepbest    number of sweep in which fbest was updated for the last
%		time
% optlevel    	level(s) of the box(es) containing the (a) global
%             	minimum of a test function
% record(1:smax-1) record(i) points to the best non-split box at level i
%             	(record list)
% xglob(1:n,1:nglob)  xglob(:,i), i=1:nglob, are the global minimizers
% of a test function in [u,v]
% xloc(1:n,:)	(for local ~= 0) columns are the points that have been 
%		used as starting points for local search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(u);

% check box bounds
if ~isempty(find(v<u)) 
  error('incompatible box bounds')
elseif ~isempty(find(u==v))
  error('degenerate box bound')
end


% default values for the input parameters 
if nargin < 5, prt = 1; end
if nargin < 6, smax = 5*n+10; end
if nargin < 7, nf = 50*n^2; end
if nargin < 8, stop = 3*n; end
if nargin < 9, 
  if isempty(find(isinf(u))) & isempty(find(isinf(v)))
    iinit = 0; 
  else
    iinit = 1;
  end
end
if nargin < 10, local = 50; end
if nargin < 11, gamma = eps; end
if nargin < 12, hess = ones(n,n); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin=[];fmi=[];		% avoid warnings in Matlab5

% initial values for the numbers of function calls (total number/local 
% search)
ncall = 0;
ncloc = 0;

% some parameters needed for initializing large arrays
step1 = 10000;
step = 1000;
dim = step1;

% initialization of some large arrays
isplit = zeros(1,step1);
level = zeros(1,step1);
ipar = zeros(1,step1);
ichild = zeros(1,step1);
f = zeros(2,step1);
z = zeros(2,step1);
nogain = zeros(1,step1);

% definition of the initialization list
if isstruct(iinit)
  x0 = iinit.x0;
  l = iinit.l;
  L = iinit.L;
  iinit = Inf;
elseif iinit == 0 
  x0(:,1) = u;
  x0(:,2) = (u+v)/2;
  x0(:,3) = v;
  l = 2*ones(n,1);
  L = 3*ones(n,1);
elseif iinit == 1
  for i = 1:n
    if u(i) >= 0
      x0(i,1) = u(i); [x0(i,2),x0(i,3)] = subint(u(i),v(i));x0(i,2) = 0.5*(x0(i,1)+x0(i,3));
    elseif v(i) <= 0
      x0(i,3) = v(i); [x0(i,2),x0(i,1)] = subint(v(i),u(i));x0(i,2) = 0.5*(x0(i,1)+x0(i,3));
    else
      x0(i,2) = 0; [xi,x0(i,1)] = subint(0,u(i)); [xi,x0(i,3)] = subint(0,v(i));
    end
  end
  l = 2*ones(n,1);
  L = 3*ones(n,1);
elseif iinit == 2
  x0(:,1) = (5*u + v)/6;
  x0(:,2) = 0.5*(u + v);
  x0(:,3) = (u + 5*v)/6;
  l = 2*ones(n,1);
  L = 3*ones(n,1);
elseif iinit == 3
  [x0,f0,l,L,istar,ncall1] = initlist(fcn,data,u,v);
  ncall = ncall + ncall1;
else
  init0 	%self-defined initialization list
  for i=1:size(x0,2)
    if ~isempty(find(x0(:,i)<u)) | ~isempty(find(x0(:,i)>v))
      error('incorrect initialization list')
    end
  end
end 

% check whether there are infinities in the initialization list
if ~isempty(find(isinf(x0))), error('infinities in ititialization list'), end

% computation of the function values f0 appertaining to the init. list 
% and the pointer istar to the best point in the initialization list
if iinit ~= 3
  [f0,istar,ncall1] = init(fcn,data,x0,l,L,n);
  ncall = ncall + ncall1; 
end

% definition of the base vertex of the original box
for i = 1:n
  x(i) = x0(i,l(i));
end

% definition of the opposite vertex v1 of the original box
for i = 1:n
  if abs(x(i)-u(i)) > abs(x(i)-v(i))
    v1(i) = u(i);
  else
    v1(i) = v(i);
  end
end

% initialization of the record list, the counters nboxes, nbasket, m 
% and nloc, xloc and the output flag
record = zeros(smax-1,1);
nboxes = 1;
nbasket = 0;
nbasket0 = 0;
nsweep = 0;
m = n;
record(1) = 1;
nloc = 0;
xloc = [];
flag = 1; 

[ipar,level,ichild,f,isplit,p,xbest,fbest] = initbox(x0,f0,l,L,istar,u,v,prt);
% generates the boxes in the initialization procedure
f0min = fbest;
if stop(1) > 0 & stop(1) < 1
  flag = chrelerr(fbest,stop);
elseif stop(1) == 0
  flag = chvtr(fbest,stop(2));
end
if ~flag,return,end
% if the (known) minimum function value fglob has been found with the
% required tolerance, flag is set to 0 and the program is terminated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = strtsw(smax,level,f(1,:)); 
% the vector record is updated, and the minimal level s containing 
% non-split boxes is computed
nsweep = nsweep + 1;	% sweep counter
  
while s < smax & ncall + 1 <= nf
  par = record(s);   % the best box at level s is the current box
  [n0,x,y,x1,x2,f1,f2] = vertex(par,n,u,v,v1,x0,f0,ipar,isplit,ichild,z,f,l,L); 
  % compute the base vertex x, the opposite vertex y, the 'neighboring' 
  % vertices and their function values needed for quadratic 
  % interpolation and the vector n0 indicating that the ith coordinate
  % has been split n0(i) times in the history of the box
  if s > 2*n*(min(n0)+1) 
  % s 'large' 
    [isplit(par),z(2,par)] = splrnk(n,n0,p,x,y);  
    % splitting index and splitting value z(2,par) for splitting by 
    % rank are computed
    % z(2,par) is set to Inf if we split according to the init. list
    splt = 1;  % indicates that the box is to be split
  else
    if nogain(par) % box has already been marked as not eligible for splitting
                   % by expected gain
      splt = 0;
    else
      [e,isplit(par),z(2,par)] = exgain(n,n0,l,L,x,y,x1,x2,f(1,par),f0,f1,f2);
      % splitting by expected gain
      % compute the expected gain vector e and the potential splitting 
      % index and splitting value
      fexp = f(1,par) + min(e);
      if fexp < fbest 
        splt = 1;
      else
        splt = 0;  % the box is not split since we expect no improvement
        nogain(par) = 1; % the box is marked as not eligible for splitting by expected gain
      end
    end
  end
  if splt == 1  % prepare for splitting
    i = isplit(par);
    level(par) = 0;
    if z(2,par) == Inf % prepare for splitting by initialization list
      m = m + 1;
      z(2,par) = m; 
      [xbest,fbest,f0(:,m),xmin,fmi,ipar,level,ichild,f,flag,ncall1] = splinit(fcn,data,i,s,smax,par,x0,n0,u,v,x,y,x1,x2,L,l,xmin,fmi,ipar,level,ichild,f,xbest,fbest,stop,prt);
      ncall = ncall + ncall1;
    else  % prepare for default splitting
      z(1,par) = x(i);
      [xbest,fbest,xmin,fmi,ipar,level,ichild,f,flag,ncall1] = split(fcn,data,i,s,smax,par,n0,u,v,x,y,x1,x2,z(:,par),xmin,fmi,ipar,level,ichild,f,xbest,fbest,stop,prt);
      ncall = ncall + ncall1;
    end
    if nboxes > dim 
% if the pre-assigned size of the `large' arrays has already been exceeded, these arrays are made larger
      isplit(nboxes+1:nboxes+step) = zeros(1,step);
      level(nboxes+1:nboxes+step) = zeros(1,step);
      ipar(nboxes+1:nboxes+step) = zeros(1,step);
      ichild(nboxes+1:nboxes+step) = zeros(1,step);
      z(:,nboxes+1:nboxes+step) = zeros(2,step);
      nogain(nboxes+1:nboxes+step) = zeros(1,step);
      f(:,nboxes+1:nboxes+step) = zeros(2,step);
      dim = nboxes + step;
    end
    if ~flag,break,end
  else  % splt=0: no splitting, increase the level by 1
    if s + 1 < smax 
      level(par) = s + 1;
      updtrec(par,s+1,f(1,:));
    else
      level(par) = 0;
      nbasket = nbasket + 1;
      xmin(:,nbasket) = x;
      fmi(nbasket) = f(1,par);
    end
    if prt > 1
      [w1,w2] = bounds(n,n0,x,y,u,v);
      % compute lower and upper bounds of the box in order to be able
      % to check whether it contains a global minimizer
      iopt = [];
      % the vector iopt contains the indices of the global minimizers
      % contained in the box
      for iglob = 1:nglob
        if w1 <= xglob(:,iglob) & xglob(:,iglob) <= w2
          iopt = [iopt, iglob];
        end
        for iglob = 1:length(iopt)
          optlevel(iopt(iglob)) = s + 1;
        end
      end      
    end
  end % of prepare for splitting
  s = s + 1;    
  while s < smax 
    if record(s) == 0
      s = s + 1;
    else
      break   
    end
  end
  if s == smax  % if smax is reached, a new sweep is started 
    if local,
      [fmi(nbasket0+1:nbasket),j] = sort(fmi(nbasket0+1:nbasket));
      xmin(:,nbasket0+1:nbasket) = xmin(:,nbasket0+j);
      xmin0 = [];
      fmi0 = [];
      for j = nbasket0+1:nbasket
        x = xmin(:,j);
        f1 = fmi(j);
        chkloc;
        if loc,
          addloc;          
          [xbest,fbest,xmin,fmi,x,f1,loc,flag,ncall1] = basket(fcn,data,x,f1,xmin,fmi,xbest,fbest,stop,nbasket0);
          ncall = ncall + ncall1;
          if ~flag,break,end
          if loc,
            [xmin1,fmi1,nc,flag] = lsearch(fcn,data,x,f1,f0min,u,v,nf-ncall,stop,local,gamma,hess);
            ncall = ncall + nc;
            ncloc = ncloc + nc;
            if fmi1 < fbest
              xbest = xmin1;
              fbest = fmi1;
              nsweepbest = nsweep;
              if ~flag
                nbasket0 = nbasket0 + 1;
                nbasket = nbasket0;
                xmin(:,nbasket) = xmin1;
                fmi(nbasket) = fmi1;
                break
              end
              if stop(1) > 0 & stop(1) < 1
                flag = chrelerr(fbest,stop);
              elseif stop(1) == 0
                flag = chvtr(fbest,stop(2));
              end
              if ~flag,return,end
            end
            [xbest,fbest,xmin,fmi,loc,flag,ncall1] = basket1(fcn,data,xmin1,fmi1,xmin,fmi,xbest,fbest,stop,nbasket0);
            ncall = ncall + ncall1;
            if ~flag,break,end
            if loc,
              nbasket0 = nbasket0 + 1;
              xmin(:,nbasket0) = xmin1;
              fmi(nbasket0) = fmi1;
              fbestloc;
              if ~flag,
                nbasket = nbasket0; break
              end
            end
          end
        end
      end
      nbasket = nbasket0;      
      if ~flag,break,end
    end
    s = strtsw(smax,level,f(1,:));
    if prt,
      if nsweep == 1
        fprintf('nsw  minl  ');
        if prt > 1
          fprintf('optl    fopt       ')
        end
        fprintf('nf     fbest        xbest\n')
      end
      minlevel=s;
      fprintf('%3i  %3i',nsweep,minlevel);
      if prt > 1
        fprintf('  %3i',optlevel);fprintf('  %10.3e',foptbox);
      end
      fprintf('  %5i  %10.3e',ncall,fbest);
      fprintf('  %10.4f',xbest);
      fprintf(1,'\n');
    end
    if stop(1) > 1
      if nsweep - nsweepbest >= stop(1),flag = 3; return,end
    end
    nsweep = nsweep + 1;
  end
end
if ncall >= nf
  flag = 2;
end
if local,
  if length(fmi) > nbasket
    xmin(:,nbasket+1:length(fmi)) = [];
    fmi(nbasket+1:length(fmi)) = [];
  end
end


