% lsquart; 	
% quartic interpolation step
% 	unfortunately, this may be unstable for huge extrapolation
%	(interpolant loses inflection points)
% 	maybe expansion around boundary points work better
%

% find quartic interpolant
if alist(1)==alist(2), f12=0;
else f12=(flist(2)-flist(1))/(alist(2)-alist(1)); 
end;
if alist(2)==alist(3), f23=0;
else f23=(flist(3)-flist(2))/(alist(3)-alist(2)); 
end;
if alist(3)==alist(4), f34=0;
else f34=(flist(4)-flist(3))/(alist(4)-alist(3)); 
end;
if alist(4)==alist(5), f45=0;
else f45=(flist(5)-flist(4))/(alist(5)-alist(4)); 
end;
f123=(f23-f12)/(alist(3)-alist(1)); 
f234=(f34-f23)/(alist(4)-alist(2));
f345=(f45-f34)/(alist(5)-alist(3));
f1234=(f234-f123)/(alist(4)-alist(1));
f2345=(f345-f234)/(alist(5)-alist(2));
f12345=(f2345-f1234)/(alist(5)-alist(1));


if f12345<=0, 
  % quartic not bounded below
  if prt>1, disp('local step (quartic not bounded below)'); end;
  good=0;lslocal;
  quart=0;
else
  if prt>1, disp('quartic step'); end; 
  quart=1;
end;

if quart,
  % expand around alist(3)
  c(1)=f12345;
  c(2)=f1234+c(1)*(alist(3)-alist(1));
  c(3)=f234+c(2)*(alist(3)-alist(4));
  c(2)=c(2)+c(1)*(alist(3)-alist(4));
  c(4)=f23+c(3)*(alist(3)-alist(2));
  c(3)=c(3)+c(2)*(alist(3)-alist(2));
  c(2)=c(2)+c(1)*(alist(3)-alist(2));
  c(5)=flist(3);
  if prt>3, 
    test_quartic_fit=[flist;quartic(c,alist-alist(3))]
  end;
  % find critical points of quartic as zeros of gradient
  cmax=max(c);c=c/cmax;
  hk=4*c(1);
  compmat= [ 0 0 -c(4); hk 0 -2*c(3); 0 hk -3*c(2)];
  ev=eig(compmat)/hk;
  i=find(imag(ev)==0);

  if max(size(i))==1,
    % only one minimizer
    if prt>1, disp('quartic has only one minimizer'); end;
    alp=alist(3)+ev(i);
    if prt>3,
      f=quartic(c,ev(i));
      plot(alp,f,'ro');
    end;
  else
    % two minimizers
    ev=sort(ev);
    if prt>1, disp('quartic has two minimizers'); end;
    alp1=lsguard(alist(3)+ev(1),alist,amax,amin,small);
    alp2=lsguard(alist(3)+ev(3),alist,amax,amin,small);
    f1=cmax*quartic(c,alp1-alist(3));
    f2=cmax*quartic(c,alp2-alist(3));
    if prt>3,
      alp3=alist(3)+ev(2);
      f3=cmax*quartic(c,ev(2));
      plot(alp1,f1,'ro');
      plot(alp2,f2,'ro');
      plot(alp3,f3,'ro');
    end;
    % pick extrapolating minimizer if possible 
    if alp2>alist(5) & f2<max(flist); alp=alp2;
    elseif alp1<alist(1) & f1<max(flist); alp=alp1;
    elseif f2<=f1, alp=alp2;
    else alp=alp1;
    end;
  end;
   
  if max(alist==alp);
    % predicted point already known
    if prt, disp('quartic predictor already known'); end;
    quart=0;
  end;
end;

if quart,
  alp=lsguard(alp,alist,amax,amin,small);
  % new function value
  falp=feval(func,data,x+alp*p);
  alist=[alist,alp];flist=[flist,falp];
  lssort;
end;


