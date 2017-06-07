%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function x1 = split2(x,y)
% determines a value x1 for splitting the interval [min(x,y),max(x,y)]
% is modeled on the function subint with safeguards for infinite y

function x1 = split2(x,y)
x2 = y;
if x == 0 & abs(y) > 1000
 x2 = sign(y);
elseif x ~= 0 & abs(y) > 100*abs(x)
 x2 = 10.*sign(y)*abs(x);
end
x1 = x + 2*(x2 - x)/3;
