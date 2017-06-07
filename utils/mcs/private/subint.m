%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subint.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2] = subint(x,y)
% computes for real x and real or infinite y two points x1 and x2 in
% [min(x,y),max(x,y)] that are neither too close nor too far away from x

function [x1,x2] = subint(x,y)
x2 = y;
f = 1000;
if f*abs(x) < 1
 if abs(y) > f, x2 = sign(y); end
else
 if abs(y) > f*abs(x), x2 = 10*sign(y)*abs(x); end
end
x1 = x + (x2 - x)/10;
