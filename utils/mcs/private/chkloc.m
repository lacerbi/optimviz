%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chkloc.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a point has already been used as starting point for a
% local search; in that case loc is set to 0, otherwise, it is set to 1

loc = 1;
for k = 1:nloc
  if x == xloc(:,k)
    loc = 0;break;
  end
end
