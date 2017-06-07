% lsnew; 	
% find one new point by extrapolation or split of wide interval
%


if prt>1,
  disp(['extrapolation or split, s = ',num2str(s)]);
elseif prt>1,
  disp(['extrapolation only, s = ',num2str(s)]);
end;

if alist(1)<=amin, 
  % leftmost point already at boundary
  leftok=0;		
elseif flist(1)>=max(fmed,flist(2)), 
  % bad leftmost point
  % extrapolate only if no scale known or global search	
  leftok= ( sinit==1 | nloc>1 );		
else 
  % good interior leftmost point
  leftok=1;				
end;
if alist(s)>=amax, 
  % rightmost point already at boundary
  rightok=0;
elseif flist(s)>=max(fmed,flist(s-1)), 
  % bad rightmost point
  % extrapolate only if no scale known or global search	
  rightok= ( sinit==1 | nloc>1 );		
else 
  % good interior rightmost point
  rightok=1;	
end;

% number of intervals used in extrapolation step
if sinit==1, step=s-1; else step=1; end;

% do the step
if leftok & ( flist(1)<flist(s) | ~rightok ),
  if prt>1, disp('extrapolate at left end point'); end;
  extra=1;
  al=alist(1)-(alist(1+step)-alist(1))/small;
  alp=max(amin,al);
elseif rightok,
  if prt>1, disp('extrapolate at right end point'); end;
  extra=1;
  au=alist(s)+(alist(s)-alist(s-step))/small;
  alp=min(au,amax);
else,
  % no extrapolation
  if prt>1, disp('split relatively widest interval'); end;
  extra=0;
  len=(alist(2:s)-alist(1:s-1));
  dist=max([alist(2:s)-abest;abest-alist(1:s-1);unitlen*ones(1,s-1)]);
  [wid,i]=max(len./dist);
  lssplit;
end;

% new function value
falp=feval(func,data,x+alp*p);
alist=[alist,alp];flist=[flist,falp];
lssort;

