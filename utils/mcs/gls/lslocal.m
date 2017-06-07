% lslocal; 		
% refine nloc best local minimizers
%


% find all strict local minima
fmin=min(flist);
up=[flist(1:s-1)<flist(2:s)];
down=[flist(2:s)<=flist(1:s-1)];
down(s-1)=(flist(s)<flist(s-1));
minima= ( [up,1] & [1,down] );
imin=find(minima);

% consider nloc best local minima only
[ff,perm]=sort(flist(imin));
imin=imin(perm);
nind=min(nloc,size(imin,2));
imin=imin(nind:-1:1);	% best point last for final improvement
if prt>3, 
  disp([alist;flist]); 
  disp([num2str(nloc),' best local minimizers at [a;f]:']);
  disp([imin;alist(imin);flist(imin)]);
elseif prt>2,
  disp([num2str(nloc),' best local minimizers at [a;f]:']);
  disp([alist(imin);flist(imin)]);
end;

nadd=0;			% number of added points
nsat=0;			% number of saturated points

for i=imin,

  % select nearest five points for local formula
  if i<=2, 
    ind=[1:5];ii=i;
  elseif i>=s-1, 
    ind=[s-4:s];ii=i-s+5;
  else 
    ind=[i-2:i+2];ii=3;
  end;
  aa=alist(ind);ff=flist(ind); 	% in natural order
  				% the local minimum is at ff(ii)

  % get divided differences 
  f12=(ff(2)-ff(1))/(aa(2)-aa(1));
  f23=(ff(3)-ff(2))/(aa(3)-aa(2));
  f34=(ff(4)-ff(3))/(aa(4)-aa(3));
  f45=(ff(5)-ff(4))/(aa(5)-aa(4));
  f123=(f23-f12)/(aa(3)-aa(1));
  f234=(f34-f23)/(aa(4)-aa(2));
  f345=(f45-f34)/(aa(5)-aa(3));
    
  % decide on action
  % cas=-1: 	no local refinement at boundary
  % cas=0: 	use parabolic minimizer
  % cas=1: 	use higher order predictor in i-2:i+1
  % cas=5: 	use higher order predictor in i-1:i+2
  % select formula on convex interval
  if ii==1,		% boundary minimum
    % parabolic minimizer or extrapolation step
    cas=0;
    if f123>0 & f123<inf, 
      alp=0.5*(aa(2)+aa(3)-f23/f123);
      if alp<amin, cas=-1; end;
    else
      alp=-inf;
      if alist(1)==amin & flist(2)<flist(3), cas=-1; end;
    end;
    alp=lsguard(alp,alist,amax,amin,small);
  elseif ii==5,		% boundary minimum
    % parabolic minimizer or extrapolation step
    cas=0;
    if f345>0 & f345<inf, 
      alp=0.5*(aa(3)+aa(4)-f34/f345);
      if alp>amax, cas=-1; end;
    else 
      alp=+inf;
      if alist(s)==amax & flist(s-1)<flist(s-2), cas=-1; end;
    end;
    alp=lsguard(alp,alist,amax,amin,small);
  elseif ~(f234>0 & f234<inf),
    % parabolic minimizer
    cas=0;
    if ii<3,
      alp=0.5*(aa(2)+aa(3)-f23/f123);
    else
      alp=0.5*(aa(3)+aa(4)-f34/f345);
    end;
  elseif ~(f123>0 & f123<inf),
    if f345>0 & f345<inf,
      cas=5; 		% use 2345
    else
      % parabolic minimizer
      cas=0;
      alp=0.5*(aa(3)+aa(4)-f34/f234);
    end;  
  elseif f345>0 & f345<inf & ff(2)>ff(4),
    cas=5; 		% use 2345
  else
    cas=1;		% use 1234
  end;

  if cas==0,
    % parabolic minimizer might extrapolate at the boundary
    alp=max(amin,min(alp,amax));
  elseif cas==1,		
    % higher order minimizer using 1234
    if ff(2)<ff(3), 
      % compute f1x4=f134
      f13=(ff(3)-ff(1))/(aa(3)-aa(1));
      f1x4=(f34-f13)/(aa(4)-aa(1));
    else 
      % compute f1x4=f124
      f24=(ff(4)-ff(2))/(aa(4)-aa(2));
      f1x4=(f24-f12)/(aa(4)-aa(1));
    end;
    alp=0.5*(aa(2)+aa(3)-f23/(f123+f234-f1x4));
    if alp<=min(aa) | alp>= max(aa),
      cas=0;
      alp=0.5*(aa(2)+aa(3)-f23/max(f123,f234));
      if prt>1, disp('predictor outside interval'); end;
    end;
  elseif cas==5, 	
    % higher order minimizer using 2345
    if ff(3)<ff(4), 
      % compute f2x5=f245
      f24=(ff(4)-ff(2))/(aa(4)-aa(2));
      f2x5=(f45-f24)/(aa(5)-aa(2));
    else 
      % compute f2x5=f235
      f35=(ff(5)-ff(3))/(aa(5)-aa(3));
      f2x5=(f35-f23)/(aa(5)-aa(2));
    end;
    alp=0.5*(aa(3)+aa(4)-f34/(f234+f345-f2x5));
    if alp<=min(aa) | alp>= max(aa),
      cas=0;
      alp=0.5*(aa(3)+aa(4)-f34/max(f234,f345));
    end;
  end;

 
  if prt>2,
    if cas>0,
      disp('higher order step');
    elseif cas==0,
      disp('parabolic step');
    else
      disp('no local refinement at boundary');
    end; 
  end;
 
  % tolerance for accepting new step
  if cas<0 | flist(i)>fmed,
    alptol=0;
  elseif cas>=0 ,
    if i==1, alptol=small*(alist(3)-alist(1));
    elseif i==s, alptol=small*(alist(s)-alist(s-2));
    else alptol=small*(alist(i+1)-alist(i-1));
    end;
  end;
  close= ( min(abs(alist-alp))<=alptol );

  if cas<0 | close,
    nsat=nsat+1;
    if prt>2, 
      if cas<0, disp('no local refinement at boundary'); 
      elseif alptol>0, disp('predicted point close to known point'); 
      else disp('predicted point matches known point'); 
      end;
    end;
  end;
  saturated=(nsat==nind);
  % check if saturated and best point changes
  final= saturated & ~max(alist==alp);
  if cas>=0 & ( final | ~close ),
    if prt>1, disp(['add local point at alp=',num2str(alp)]); end;
    % add point to the list
    nadd=nadd+1;
    % new function value
    falp=feval(func,data,x+alp*p);
    alist=[alist,alp];flist=[flist,falp];
    % no sort since this would destroy old index set!!!
    if prt>1, abest_anew_xnew=[alist(i),alp,(x+alp*p)']; end;
  end;

end;

if nadd, lssort; end;

if prt>1,
  if saturated, disp(['saturated at s = ',num2str(s)]); 
  else disp(['not saturated at s = ',num2str(s)]); 
  end; 
end; 

