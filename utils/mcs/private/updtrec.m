%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updtrec.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function updtrec(j,s,f)
% Input:
% j           label of a box
% s           its level
% f           vector containing the base vertex function values of the
%             already defined boxes
% updates the pointer record(s) to the best non-split box at level s
function updtrec(j,s,f)
global record      % record list
if length(record) < s
  record(s) = j;
elseif record(s) == 0
  record(s) = j;
elseif f(j) < f(record(s))
  record(s) = j;
end
