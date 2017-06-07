% lssep; 		
% separate close local minimizers
% and maybe add a few global points
%


nsep=0;

while nsep<nmin,
  % find intervals where the monotonicity behavior of both
  % adjacent intervals is opposite
  down=[flist(2:s)<flist(1:s-1)];
  sep= ( [1,1,down] & [0,up,0] & [down,1,1] );
  sep= sep | ( [1,1,up] & [0,down,0] & [up,1,1] );
  ind=find(sep);
  if isempty(ind), break; end; 

  aa=0.5*(alist(ind)+alist(ind-1));	% interval midpoints
  if length(aa)>nloc,
    % select nloc best interval midpoints
    ff=min(flist(ind),flist(ind-1));
    [ff,ind]=sort(ff);
    aa=aa(ind(1:nloc));
  end;
  for alp=aa,
    if prt>2, disp(['separate minimizer at ',num2str(alp)]); end;
    % new function value
    falp=feval(func,data,x+alp*p);
    alist=[alist,alp];flist=[flist,falp];
    nsep=nsep+1;
    if nsep>=nmin, break; end; 
  end;
  lssort;
end;

% instead of unnecessary separation, add some global points
for times=1:nmin-nsep, 
  lsnew;			% extrapolation or split
end;
