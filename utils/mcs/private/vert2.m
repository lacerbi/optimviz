%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% vert2.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x1,x2,f1,f2] = vert2(j,x,z,f,x1,x2,f1,f2)
% subprogram used by vertex.m

function [x1,x2,f1,f2] = vert2(j,x,z,f,x1,x2,f1,f2)
if j == 1
  j1 = 2;
else
  j1 = 1;
end
if x1 == Inf
  x1 = z(j);
  f1 = f1 + f(j);
  if x ~= z(j1)
    x2 = z(j1);
    f2 = f2 + f(j1);
  end
elseif x2 == Inf & x1 ~= z(j)
  x2 = z(j);
  f2 = f2 + f(j);
elseif x2 == Inf
  x2 = z(j1);
  f2 = f2 + f(j1);
end
