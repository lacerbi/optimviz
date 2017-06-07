%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chrelerr.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function flag = chrelerr(fbest,stop)
% checks whether the required tolerance for a test function with known
% global minimum has already been achieved
% Input:
% fbest		function value to be checked
% stop(1)	relative error with which a global minimum with not too
%		small absolute value should be reached
% stop(2)	global minimum function value of a test function
% stop(3)	if abs(fglob) is very small, we stop if the function
%		value is less than stop(3)
% Output:
% flag          = 0 the required tolerance has been achieved
% 		= 1 otherwise

function flag = chrelerr(fbest,stop)
fglob = stop(2);
if fbest - fglob <= max(stop(1)*abs(fglob),stop(3))
  flag = 0;
else
  flag = 1;
end   

