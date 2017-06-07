%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hessian.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function h = mcs_hessian(i,k,x,x0,f,f0,g,G)
% computes the element G(i,k) of the Hessian of the local quadratic 
% model
% Input:
% i, k		indices (k < i)
% x(1:n)	'neighbor' of x0 used for computing G(i,k); differs from
%		x0 only in the ith and kth component
% x0(1:n)	point around which the quadratic model is computed
% f		function value at x
% f0		function value at x0
% g(1:i)	components of the gradient of the local quadratic model %		that have already been computed
% G		components of the Hessian of the local quadratic model
%		that have already been computed
% Output:
% h = G(i,k)	newly computed nondiagonal element of the Hessian

function h = mcs_hessian(i,k,x,x0,f,f0,g,G)
h = f-f0-g(i)*(x(i)-x0(i))-g(k)*(x(k)-x0(k))-0.5*G(i,i)*(x(i)-x0(i))^2-0.5*G(k,k)*(x(k)-x0(k))^2;
h = h/(x(i)-x0(i))/(x(k)-x0(k));
