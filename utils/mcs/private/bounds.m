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
function [w1,w2] = bounds(n,n0,x,y,u,v)
w1 = zeros(n,1);
w2 = w1;
for i = 1:n
  if n0(i) == 0
    w1(i) = u(i);
    w2(i) = v(i);
  else
    w1(i) = min(x(i),y(i));
    w2(i) = max(x(i),y(i));
  end
end
