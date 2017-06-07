% lspar; 		
% local parabolic minimization
%

cont=1;	% continue?
if s<3, lsnew; cont=0; end;

if cont,
  % select three points for parabolic interpolation
  [fmin,i]=min(flist);
  if i<=2 , 
    ind=[1:3];ii=i;
  elseif i>=s-1, 
    ind=[s-2:s];ii=i-s+3;
  else 
    ind=[i-1:i+1];ii=2;
  end;
  aa=alist(ind);ff=flist(ind); 	% in natural order
				% the local minimum is at ff(ii)

  % get divided differences 
  f12=(ff(2)-ff(1))/(aa(2)-aa(1));
  f23=(ff(3)-ff(2))/(aa(3)-aa(2));
  f123=(f23-f12)/(aa(3)-aa(1));		

  % handle concave case
  if ~(f123>0),
    if prt>1,disp('parabola not convex');end; 
    if prt>3,f123,end; 
    lsnew;cont=0; 
  end;
end;

if cont,
  % parabolic minimizer
  alp0=0.5*(aa(2)+aa(3)-f23/f123);
  alp=lsguard(alp0,alist,amax,amin,small);
  alptol=small*(aa(3)-aa(1));

  % handle infinities and close predictor
  if f123==inf | min(abs(alist-alp))<=alptol, 
    % split best interval
    if prt>1, disp('split best interval'); end;
    if ii==1 | ( ii==2 & aa(2)>=0.5*(aa(1)+aa(3)) ),
      alp=0.5*(aa(1)+aa(2));
    else 
      alp=0.5*(aa(2)+aa(3));
    end;
  else
    if prt>1, disp('parabolic step'); end;
    if prt>3 & alp~=alp0, 
      disp(['parabolic predictor: alp0 = ',num2str(alp0)]); 
    end; 
  end;

  % add point to the list
  if prt>1, disp(['add point at alp=',num2str(alp)]); end;
  % new function value
  falp=feval(func,data,x+alp*p);
  alist=[alist,alp];flist=[flist,falp];
  if prt>1, abest_anew_xnew=[alist(i),alp,(x+alp*p)']; end;
end;

lssort;


