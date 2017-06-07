%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtoptl.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function updtoptl(i,x,y,iopt,level,f)
% a daughter box of a box containing (a) local minimizer(s) is checked
% whether it contains any minimizers, and in case it contains any, the
% corresponding levels and function values are updated 
% Input:
% i      splitting index of the parent box
% x,y    bounds of the ith coordinate of the daughter box 
% iopt   vector containing the indices of the global minimizers the
%        parent box contained
% level  level of the daughter box
% f      base function value of the daughter box
function updtoptl(i,x,y,iopt,level,f)
global foptbox optlevel xglob
% foptbox(1:nglob)   base vertex function value(s) of the box(es) 
%        containing the (a) global minimizer of a test function
% optlevel(1:nglob)  their level(s)
% xglob(1:n,1:nglob) xglob(:,i), i=1:nglob, are the global minimizers 
% of a test function 
for j = 1:length(iopt)
  if min(x,y) <= xglob(i,iopt(j)) & xglob(i,iopt(j)) <= max(x,y)
    optlevel(iopt(j)) = level;
    foptbox(iopt(j)) = f; 
  end
end
