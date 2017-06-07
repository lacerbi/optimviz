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

function [g,G] = polint1(x,f)
f13 = (f(3)-f(1))/(x(3)-x(1));
f12 = (f(2)-f(1))/(x(2)-x(1));
f23 = (f(3)-f(2))/(x(3)-x(2));
g = f13 + f12 - f23;
G = 2*(f13-f12)/(x(3)-x(2));
