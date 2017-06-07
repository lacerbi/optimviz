% lssat;		
% check saturation condition
% 

cont=saturated;

if cont, 
  % check boundary minimizer
  [fmin,i]=min(flist);
  if i==1 | i==s, cont=0; end;
end;

if cont, 
  % select three points for parabolic interpolation
  aa=alist([i-1:i+1]);ff=flist([i-1:i+1]); 	

  % get divided differences 
  f12=(ff(2)-ff(1))/(aa(2)-aa(1));
  f23=(ff(3)-ff(2))/(aa(3)-aa(2));
  f123=(f23-f12)/(aa(3)-aa(1));

  if f123>0,
    % parabolic minimizer
    alp=0.5*(aa(2)+aa(3)-f23/f123);
    alp=max(amin,min(alp,amax));
    alptol=small*(aa(3)-aa(1));
    saturated= ( abs(alist(i)-alp)<=alptol );
  else
    saturated=0;
  end;
  
  if prt>1 & ~saturated,
    disp('saturation check negative')
  end;
end;
