% lssort;		
% sort points by increasing alp;
% count list size and number of strict local extrema
%
% simple but inefficient 
% may be updated in a much faster way 
% by looking at in which intervals points are added
%

[alist,perm]=sort(alist);flist=flist(perm);
s=size(alist,2);

% find number of strict local minima, etc.
up=[flist(1:s-1)<flist(2:s)];
down=[flist(2:s)<=flist(1:s-1)];
down(s-1)=(flist(s)<flist(s-1));
monotone= ( sum(up)==0 | sum(down)==0 );
minima=([up,1] & [1,down]);
nmin=sum(minima);
[fbest,i]=min(flist);abest=alist(i);
fmed=median(flist);
 
% distance from best minimum to next
if nmin>1,
  al=alist(minima);al(al==abest)=[];
  unitlen=min(abs(al-abest));
else
  unitlen=max(abest-alist(1),alist(s)-abest);
end;

