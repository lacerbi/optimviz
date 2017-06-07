

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% lsguard.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function alp=lsguard(alp,alist,amax,amin,small); 		
% safeguard a point alp
% alist			current list of steps
% amax,amin		admissible range for alp
% small			relative minimal distance from end point
% 
% assumes that alist(1:s) is within these bounds

function alp=lsguard(alp,alist,amax,amin,small); 	

asort=sort(alist);
s=size(asort,2);

% enforce extrapolation to be cautious
al=asort(1)-(asort(s)-asort(1))/small;
au=asort(s)+(asort(s)-asort(1))/small;
alp=max(al,min(alp,au));
alp=max(amin,min(alp,amax));

% enforce some distance from end points
% 	factor 1/3 ensures equal spacing if s=2 and the third point
% 	in a safeguarded extrapolation is the maximum.
if abs(alp-asort(1))<small*(asort(2)-asort(1)),
  alp=(2*asort(1)+asort(2))/3;
end;
if abs(alp-asort(s))<small*(asort(s)-asort(s-1)),
  alp=(2*asort(s)+asort(s-1))/3;
end;
