% lsinit; 	
% find first two points and establish correct scale
%

if size(alist,2)==0,
  % evaluate at absolutely smallest point
  alp=0;
  if amin>0, alp=amin; end;
  if amax<0, alp=amax; end;
  % new function value
  falp=feval(func,data,x+alp*p);
  alist=alp;flist=falp;
elseif size(alist,2)==1,
  % evaluate at absolutely smallest point
  alp=0;
  if amin>0, alp=amin; end;
  if amax<0, alp=amax; end;
  if alist~=alp,
    % new function value
    falp=feval(func,data,x+alp*p);
    alist=[alist,alp];flist=[flist,falp];
  end;
end;
aamin=min(alist);aamax=max(alist);
if amin>aamin | amax<aamax,
  alist,amin,amax
  error('non-admissible step in alist');
end;

% establish correct scale
if aamax-aamin<=scale,
  alp1=max(amin,min(-scale,amax));
  alp2=max(amin,min(+scale,amax));
  alp=inf;
  if aamin-alp1>=alp2-aamax, alp=alp1; end;  
  if alp2-aamax>=aamin-alp1, alp=alp2; end;  
  if alp<aamin | alp>aamax,
    % new function value
    falp=feval(func,data,x+alp*p);
    alist=[alist,alp];flist=[flist,falp];
  end;
end;

if size(alist,2)==1,
  scale,aamin,aamax,alp1,alp2
  error('lsinit bug: no second point found'); 
end;

lssort;
