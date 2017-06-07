% function [alist,flist,gTp,pTGp,nmin,new]=lspost(alist,flist,prt) 		
% postprocessing of line search output
% find local quadratic model around best point
% and discard nonminimal points
%
% alist,flist:	list of steps and corresponding function values 
%		at input: steps with known function values
%		at output: extrema found, except bad ones at boundary
		% sorted by increasing f
%		in particular, the best point is at alist(1)
% prt:		print level
%		0: nothing plotted or printed (default)
% 		1: mark extrema in plot
% gTp,pTGp	gradient and Hessian information from 
%		local quadratic model around best point
% nmin:		number of local minima detected
% new:		is alp=0 separated from abest by a barrier?
%

%
function [alist,flist,gTp,pTGp,nmin,new]=lspost(alist,flist,prt)

if nargin<3, prt=0; end;
s=size(alist,2);

% calculate local quadratic model information
[fbest,i]=min(flist);
if i==1, ind=[1,2,3];
elseif i==s, ind=[s,s-1,s-2];
else ind=[i,i-1,i+1];
end;
aa=alist(ind);ff=flist(ind);
f12=(ff(2)-ff(1))/(aa(2)-aa(1));
f13=(ff(3)-ff(1))/(aa(3)-aa(1));
f23=(ff(3)-ff(2))/(aa(3)-aa(2));
f123=(f23-f12)/(aa(3)-aa(1));
gTp=f12+f13-f23;
pTGp=2*f123;

% find local extrema
up=[flist(1:s-1)<=flist(2:s)];		% <= needed for constant f
down=[flist(2:s)<=flist(1:s-1)];	% <= needed for constant f
maxima=([down,1] & [1,up]);
minima=([up,1] & [1,down]);
nmin=sum(minima);


if prt, 
  figure(1);
  plot(alist(maxima),flist(maxima),'+');
  plot(alist(minima),flist(minima),'+');drawnow;
  pause(0.5);
end;

% purge list from nonextremal points
ind=(maxima | minima);
if ind==[], error('1???'); end; 
alist=alist(ind);flist=flist(ind);
s=size(alist,2);

% check for separating barrier
% (sitting on a local maximum is also counted as barrier)
i=max(find(alist<=0));
if i==[], f0=flist(1);
elseif i==s, f0=flist(s);
else f0=min(flist(i),flist(i+1));
end;
new= (f0>min(flist));


% sort list by function value
[flist,perm]=sort(flist);
alist=alist(perm);

% discard bad boundary points
for k=1:2,
  if s==0, error('2???'); end; 
  if s<=1, break; end;
  if alist(s)==min(alist) | alist(s)==max(alist),
    alist(s)=[];flist(s)=[];s=s-1;
  else
    break;
  end;
end;

