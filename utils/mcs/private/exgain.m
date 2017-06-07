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
%
% Uses the following functions/m-files:
% polint.m
% quadmin.m
% quadpol.m
% subint.m

function [e,isplit,splval] = exgain(n,n0,l,L,x,y,x1,x2,fx,f0,f1,f2)
emin = Inf;  % initialization
for i = 1:n
  if n0(i) == 0
    e(i) = min(f0(1:L(i),i)) - f0(l(i),i);  
    % expected gain for splitting according to the initialization list
    if e(i) < emin
      emin = e(i);
      isplit = i;
      splval = Inf;
    end
  else
    z1 = [x(i) x1(i) x2(i)];
    z2 = [0 f1(i) - fx f2(i) - fx];
    d = polint(z1,z2);
    [eta1,eta2] = subint(x(i),y(i));
    % safeguard against splitting too close to x(i)
    xi1 = min(eta1,eta2);
    xi2 = max(eta1,eta2);
    z = quadmin(xi1,xi2,d,z1);
    e(i) = quadpol(z,d,z1);
    if e(i) < emin
      emin = e(i);
      isplit = i;
      splval = z;
    end
  end
end

