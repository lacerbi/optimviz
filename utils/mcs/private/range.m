

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% range.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function range(u,v,x,p)
% computes the range [amin,amax] for making a line search along x + ap
% when x + ap is restricted to [u,v]
%
function [amin,amax] = range(u,v,x,p)
amin = -Inf; amax = Inf;
for i=1:length(x)
  if p(i)>0, 
    amin=max(amin,(u(i)-x(i))/p(i));
    amax=min(amax,(v(i)-x(i))/p(i));
  elseif p(i)<0, 
    amin=max(amin,(v(i)-x(i))/p(i));
    amax=min(amax,(u(i)-x(i))/p(i));
  end
end  
