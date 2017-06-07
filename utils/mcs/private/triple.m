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
% Uses the following m-files:
% hessian.m
% polint1.m

function [xtrip,ftrip,g,G,x1,x2,nf] = triple(fcn,data,x,f,x1,x2,u,v,hess,G)
nf = 0;
n = length(x);
g = zeros(n,1);
if nargin < 9
  hess = ones(n,n);
end
ind = find(u < x & x < v);
ind1 = find(x <= u | x >= v);
for j=1:length(ind1)
  g(ind1(j)) = 0;
  for k=1:n
    G(ind1(j),k) = 0;
    G(k,ind1(j)) = 0;
  end
end	
if length(ind) <= 1
  xtrip = x;
  ftrip = f;
  if ~isempty(ind) 
    g(ind) = 1;
    G(ind,ind) = 1;
  end
  return
end
if nargin < 10
  G = zeros(n,n);
end	
xtrip = x;
ftrip = f;
xtripnew = x;
ftripnew = f;
for j=1:length(ind)
  i = ind(j);
  x = xtrip;
  f = ftrip;
  x(i) = x1(i);
  f1 = feval(fcn,data,x);
  x(i) = x2(i);
  f2 = feval(fcn,data,x);
  nf = nf + 2;
  [g(i),G(i,i)] = polint1([xtrip(i) x1(i) x2(i)],[f f1 f2]);
  if f1 <= f2
    if f1 < ftrip
      ftripnew = f1;
      xtripnew(i) = x1(i);
    end
  else
    if f2 < ftrip
      ftripnew = f2;
      xtripnew(i) = x2(i);
    end
  end 
  if nargin < 10
    k1 = 0;
    if f1 <= f2
      x(i) = x1(i);
    else
      x(i) = x2(i);
    end
    for k=1:i-1
      if hess(i,k)
        if xtrip(k) > u(k) & xtrip(k) < v(k) & ~isempty(find(ind==k))
          q1 = ftrip + g(k)*(x1(k)-xtrip(k))+0.5*G(k,k)*(x1(k)-xtrip(k))^2;
          q2 = ftrip + g(k)*(x2(k)-xtrip(k))+0.5*G(k,k)*(x2(k)-xtrip(k))^2; 
          if q1 <= q2
            x(k) = x1(k);
          else
            x(k) = x2(k);
          end
          f12 = feval(fcn,data,x);
          nf = nf + 1;
          G(i,k) = hessian(i,k,x,xtrip,f12,ftrip,g,G);
          G(k,i) = G(i,k);
          if f12 < ftripnew
            ftripnew = f12;
            xtripnew = x;
            k1 = k;
          end  
          x(k) = xtrip(k);
        end
      else
        G(i,k) = 0;
        G(k,i) = 0;
      end
    end
  end
  if ftripnew < ftrip
    if x1(i) == xtripnew(i)
      x1(i) = xtrip(i);
    else 
      x2(i) = xtrip(i);
    end
    if nargin < 10 & k1 > 0
      if xtripnew(k1) == x1(k1)
        x1(k1) = xtrip(k1);
      else
        x2(k1) = xtrip(k1);
      end
    end
    for k=1:i
      if ~isempty(find(ind==k))
        g(k) = g(k) + G(i,k)*(xtripnew(i) - xtrip(i));
        if nargin < 10 & k1 > 0
          g(k) = g(k) + G(k1,k)*(xtripnew(k1) - xtrip(k1));
        end
      end
    end
    xtrip = xtripnew;
    ftrip = ftripnew;
  end
end
