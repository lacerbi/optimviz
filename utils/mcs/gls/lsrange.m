% lsrange;
% find range of useful alp in truncated or bent line search
%

if max(abs(p))==0, error('zero search direction in line search'); end;

% find sensible step size scale
pp=abs(p(p~=0));
u=abs(x(p~=0))./pp; 
scale=min(u);
if scale==0, 
  u(u==0)=1./pp(u==0); 
  scale=min(u);
end;

if ~bend,
  % find range of useful alp in truncated line search
  amin=-inf;amax=+inf;
  for i=1:size(x,1),
    if p(i)>0, 
      amin=max(amin,(xl(i)-x(i))/p(i));
      amax=min(amax,(xu(i)-x(i))/p(i));
    elseif p(i)<0, 
      amin=max(amin,(xu(i)-x(i))/p(i));
      amax=min(amax,(xl(i)-x(i))/p(i));
    end;
  end;
  if amin>amax, error('no admissible step in line search'); end; 
  if prt,
    figure(1);clf;drawnow;hold on; 
% end; if 0, % this removes drawing of the curve
    aa=amin+[0:100]*(amax-amin)/100;ff=[];
    for alp=aa, 
      xx=max(xl,min(x+alp*p,xu));
      ff=[ff,feval(func,data,xx)]; 
    end;
    plot(aa,ff,'r');
    plot(alist,flist,'*');drawnow;
  end;
  if prt>1,
    disp('range of alp in truncated line search:');
    disp([amin,amax]);
  end; 
else
  % find range of useful alp in bent line search
  amin=+inf;amax=-inf;
  for i=1:size(x,1),
    if p(i)>0, 
      amin=min(amin,(xl(i)-x(i))/p(i));
      amax=max(amax,(xu(i)-x(i))/p(i));
    elseif p(i)<0, 
      amin=min(amin,(xu(i)-x(i))/p(i));
      amax=max(amax,(xl(i)-x(i))/p(i));
    end;
  end;
  if prt,
    aa=amin+[0:100]*(amax-amin)/100;ff=[];
    for alp=aa, 
      xx=max(xl,min(x+alp*p,xu));
      ff=[ff,feval(func,data,xx)]; 
    end;
    figure(1);clf; 
    plot(aa,ff,'r');hold on;
    plot(alist,flist,'*');drawnow;
  end;
  if prt>1,
    disp('range of alp in bent line search:');
    amin,amax 
  end; 
end;

if prt>2, % | (prt>0 & amin==0) | (prt>0 & amax==0), 
  disp('ls info:')
  xl,xu,x,p, 
end;
