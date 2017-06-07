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
%
% Uses the following m-files:
% chrelerr.m
% csearch.m 
% chvtr.m
% hessian.m	called by triple.m
% gls.m and its subprograms
% minq.m and its subprograms
% neighbor.m	
% polint1.m	called by triple.m
% range.m
% triple.m
function [xmin,fmi,ncall,flag] = lsearch(fcn,data,x,f,f0,u,v,nf,stop,maxstep,gamma,hess)
global nsweep nsweepbest
ncall = 0;
n = length(x);
x0 = min(max(u,0),v); % absolutely smallest point
if nargin < 10, maxstep = 50; end
if nargin < 11, gamma = eps; end
if nargin < 12, hess = ones(n,n); end
flag = 1;
eps0 = 0.001;
nloc = 1;
small = 0.1;
smaxls = 15;
[xmin,fmi,g,G,nfcsearch] = csearch(fcn,data,x,f,u,v,hess);
xmin = max(u,min(xmin,v));
ncall = ncall + nfcsearch;
xold = xmin;
fold = fmi;
if stop(1) > 0 & stop(1) < 1
  flag = chrelerr(fmi,stop);
elseif stop(1) == 0
  flag = chvtr(fmi,stop(2));
end 
if ~flag,return,end
d = min(min(xmin-u,v-xmin),0.25*(1+abs(x-x0)));
p = minq(fmi,g,G,-d,d,0);
if norm(p),
  x = xmin + p;
  x = max(u,min(x,v));
  f1 = feval(fcn,data,x);
  ncall = ncall + 1;
  alist = [0 1];
  flist = [fmi f1];
  fpred = fmi + g'*p + 0.5*p'*G*p;
  [alist,flist,nfls] = gls(fcn,data,u,v,xmin,p,alist,flist,nloc,small,smaxls);
  ncall = ncall + nfls;
  [fminew,i] = min(flist);
  if fminew == fmi
    i = find(~alist);
  else
    fmi = fminew;
  end
  xmin = xmin + alist(i)*p;
  xmin = max(u,min(xmin,v));
  gain = f - fmi;
  if stop(1) > 0 & stop(1) < 1
    flag = chrelerr(fmi,stop);
  elseif stop(1) == 0
    flag = chvtr(fmi,stop(2));
  end
  if ~flag,return,end
  if fold == fmi
    r = 0;
  elseif fold == fpred
    r = 0.5;
  else
    r = (fold-fmi)/(fold-fpred);
  end
else
  gain = f - fmi;
  r = 0;
end
diag = 0;
ind = find(u < xmin & xmin < v);
b = abs(g)'*max(abs(xmin),abs(xold));
nstep = 0;
while ncall < nf & nstep < maxstep & (diag | length(ind) < n | (stop(1) == 0 & fmi - gain <= stop(2)) | (b >= gamma*(f0-f) & gain > 0)) 
  nstep = nstep + 1;
  delta = abs(xmin)*eps^(1/3);
  j = find(~delta);
  if ~isempty(j)
    delta(j) = eps^(1/3)*ones(size(j));
  end
  [x1,x2] = neighbor(xmin,delta,u,v);
  f = fmi;
  if length(ind) < n & (b < gamma*(f0-f) | ~gain)
    ind1 = find(xmin == u | xmin == v);
    for k=1:length(ind1)
      i = ind1(k);
      x = xmin;
      if xmin(i) == u(i)
        x(i) = x2(i);
      else
        x(i) = x1(i);
      end
      f1 = feval(fcn,data,x);
      ncall = ncall + 1;
      if f1 < fmi
        alist = [0 x(i)-xmin(i)];
        flist = [fmi f1];
        p = zeros(n,1);
        p(i) = 1;
        [alist,flist,nfls] = gls(fcn,data,u,v,xmin,p,alist,flist,nloc,small,6);
        ncall = ncall + nfls;
        [fminew,j] = min(flist);
        if fminew == fmi
          j = find(~alist);
        else
          fmi = fminew;
        end
        xmin(i) = xmin(i) + alist(j);
      else
        ind1(k) = 0;
      end
    end
    xmin = max(u,min(xmin,v));
    if ~sum(ind1),break,end
    delta = abs(xmin)*eps^(1/3);
    j = find(~delta);
    if ~isempty(j)
      delta(j) = eps^(1/3)*ones(size(j));
    end
    [x1,x2] = neighbor(xmin,delta,u,v);
  end 
  if abs(r-1) > 0.25 | ~gain | b < gamma*(f0-f)
    [xmin,fmi,g,G,x1,x2,nftriple] = triple(fcn,data,xmin,fmi,x1,x2,u,v,hess);
    ncall = ncall + nftriple;
    diag = 0;
  else
    [xmin,fmi,g,G,x1,x2,nftriple] = triple(fcn,data,xmin,fmi,x1,x2,u,v,hess,G);
    ncall = ncall + nftriple;
    diag = 1;
  end
  xold = xmin;
  fold = fmi;
  if stop(1) > 0 & stop(1) < 1
    flag = chrelerr(fmi,stop);
  elseif stop(1) == 0
    flag = chvtr(fmi,stop(2));
  end
  if ~flag,return,end
  if r < 0.25
    d = 0.5*d;
  elseif r > 0.75
    d = 2*d;
  end
  p = minq(fmi,g,G,max(-d,u-xmin),min(d,v-xmin),0);
  if ~norm(p) & ~diag & length(ind) == n,break,end
  if norm(p),
    fpred = fmi + g'*p + 0.5*p'*G*p;
    x = xmin + p;
    f1 = feval(fcn,data,x);
    ncall = ncall + 1;
    alist = [0 1];
    flist = [fmi f1];
    [alist,flist,nfls] = gls(fcn,data,u,v,xmin,p,alist,flist,nloc,small,smaxls);
    ncall = ncall + nfls;
    [fmi,i] = min(flist);
    xmin = xmin + alist(i)*p;
    xmin = max(u,min(xmin,v));
    if stop(1) > 0 & stop(1) < 1
      flag = chrelerr(fmi,stop);
    elseif stop(1) == 0
      flag = chvtr(fmi,stop(2));
    end
    if ~flag,return,end 
    gain = f - fmi;
    if fold == fmi
      r = 0;
    elseif fold == fpred
      r = 0.5;
    else
      r = (fold-fmi)/(fold-fpred);
    end
    if fmi < fold  
      fac = abs(1-1/r);
      eps0 = max(eps,min(fac*eps0,0.001));
    else
      eps0 = 0.001;
    end
  else
    gain = f - fmi;
    if ~gain
      eps0 = 0.001;
      fac = Inf;
      r = 0;
    end
  end
  ind = find(u < xmin & xmin < v);
  b = abs(g)'*max(abs(xmin),abs(xold));
end

