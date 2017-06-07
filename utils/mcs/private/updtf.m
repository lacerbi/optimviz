%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtf.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [f1,f2,fold] = updtf(n,i,x1,x2,f1,f2,fold,f)
% subprogram used by vertex.m
% If we leave the i-coordinate x(i) of the base vertex of a box by
% going back to the base vertices of its ancestors, we have to 
% generate 'fictitious' function values by assuming a separable model
% for all coordinates i1 ~= i for which coordinate and function values
% have not yet been found
% Input:
% n        dimension of the problem
% i        splitting index
% x1(1:n), x2(1:n), f1(1:n), f2(1:n)
%          x1(j) and x2(j) and the corresponding function values f1(j)
%          and f2(j) used for quadratic interpolation in the jth
%          coordinate 
%          x1(j) = Inf or x2(j) = Inf indicates that these quantities
%          have not been found yet
% fold     base vertex function value of the previously considered box
% f        base vertex function value of the current box
% Output
% f1(1:n), f2(1:n) updates of f1 and f2 
% fold     base vertex function value of the current box

function [f1,f2,fold] = updtf(n,i,x1,x2,f1,f2,fold,f)
for i1=1:n
  if i1 ~= i
    if x1(i1) == Inf
      f1(i1) = f1(i1) + fold - f;
    end 
    if x2(i1) == Inf
      f2(i1) = f2(i1) + fold - f;
    end
  end
end
fold = f;
