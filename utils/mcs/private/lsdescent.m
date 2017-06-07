% lsdescent;		
% check descent condition
% 

cont=max(alist==0); % condition for continue

if cont, 
  [fbest,i]=min(flist);
  if alist(i)<0, 
    if alist(i)>=4*alist(i+1), cont=0; end; 
  elseif alist(i)>0, 
    if alist(i)<4*alist(i-1), cont=0; end; 
  else 
    if i==1, fbest=flist(2);
    elseif i==s, fbest=flist(s-1);
    else fbest=min(flist(i-1),flist(i+1));
    end;
  end;
end;

if cont, 
  % force local descent step
  if alist(i)~=0, alp=alist(i)/3;
  elseif i==s, alp=alist(s-1)/3;
  elseif i==1, alp=alist(2)/3;
  else 
    % split wider adjacent interval
    if alist(i+1)-alist(i)>alist(i)-alist(i-1), alp=alist(i+1)/3;
    else alp=alist(i-1)/3;
    end;
  end;

  % new function value
  falp=feval(func,data,x+alp*p);
  alist=[alist,alp];flist=[flist,falp];
  lssort;
  if prt>1, disp(['descent check: new point at ',num2str(alp)]); end;
end;
  
