%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% init.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f0,istar,ncall] = init(fcn,data,x0,l,L,n)
% computes the function values corresponding to the initialization list
% and the pointer istar to the final best point x^* of the init. list
% Input:
% fcn = 'fun' 	name function fun(data,x), x an n-vector
% data	      	data vector 
% x0(1:n,1:max(L))  coordinates used in the initialization list  
% l(1:n)      	pointer to the initial point x^0 = (x0(1,l(1)),...,
%             	x0(n,l(n)) of the initialization list
% L(1:n)      	the init. list contains L(i) values for coordinate i
% n           	dimension of the problem
% Output:
% f0(1:max(L),1:n)  function values appertaining to the initialization
%             	list; f0(j,i) = function value at x = (x0(1,istar(1)),
%		...,x0(i-1,istar(i-1)),x0(i,j),x0(i+1,l(i+1)),...,
%		x0(n,l(n)))
% istar(1:n)  	pointer to the final best point x^* of the init. list
%             	x^* = (x0(1,istar(1)),...,x0(n,istar(n)))
% ncall	      	number of function calls used in the program

function [f0,istar,ncall] = init(fcn,data,x0,l,L,n)
ncall = 0;
for i = 1:n
  x(i) = x0(i,l(i));
end
x = x';
f1 = feval(fcn,data,x);
f0(l(1),1) = f1;
ncall = ncall + 1; 
for i = 1:n
  istar(i) = l(i);
  for j = 1:L(i)
    if j == l(i)
      if i ~= 1
        f0(j,i) = f0(istar(i-1),i-1);
      end
    else
      x(i) = x0(i,j);
      f0(j,i) = feval(fcn,data,x);
      ncall = ncall + 1;
      if f0(j,i) < f1
        f1 = f0(j,i);
        istar(i) = j;
      end
    end
  end
  x(i) = x0(i,istar(i));
end
