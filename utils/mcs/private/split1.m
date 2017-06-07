%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% split1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function z = split1(x1,x2,f1,f2)
% Input: points x1 and x2, x1 < x2, and corresponding function values f1
%        and f2
% splits the interval [x1,x2] according to the golden section rule
% the part containing the better point gets the larger fraction of the 
% interval

function z = split1(x1,x2,f1,f2)
if f1 <= f2
  z = x1 + 0.5*(-1 + sqrt(5))*(x2 - x1);
else
  z = x1 + 0.5*(3 - sqrt(5))*(x2 - x1);
end
