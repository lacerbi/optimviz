%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fbestloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a new point in the shopping basket (found by local 
% search) is better than the current best point; in that case it 
% updates the best point and its function value
% Uses the following m-file:
% chrelerr.m

if fmi(nbasket0) < fbest
  fbest = fmi(nbasket0)
  xbest = xmin(:,nbasket0)
  chrelerr;
end
