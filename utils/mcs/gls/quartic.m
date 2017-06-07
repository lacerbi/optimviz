% function f=quartic(a,x);
% evaluates f(x)=a(1)x^4+a(2)x^3+a(3)x^2+a(4)x+a(5)
% simultaneously at a vector x
%
function f=quartic(a,x);

f=(((a(1)*x+a(2)).*x+a(3)).*x+a(4)).*x+a(5); 

