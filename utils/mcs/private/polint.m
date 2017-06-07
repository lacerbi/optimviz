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

function d = polint(x,f)
d(1) = f(1);
d(2) = (f(2) - f(1))/(x(2) - x(1));
f23 = (f(3) - f(2))/(x(3) - x(2));
d(3) = (f23 - d(2))/(x(3) - x(1));
