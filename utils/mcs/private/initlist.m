%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initlist.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x0,f0,l,L,istar,ncall] = initlist(fcn,data,u,v)
% generates an initialization list with the aid of line searches
%
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector
% [u,v]       	box in which the optimization is carried out (u, v 
%             	n-vectors)
% Output:
% x0		array with n rows and at least 3 columns; the ith row
%		contains the initialization list values for the ith 
%		coordinate 
% f0		corresponding function values
% l		x0(i,l(i)) is the ith coordinate of the initial point
% L		
% istar		x0(i,istar(i)) is the ith coordinate of the final best 
%		point
% ncall		number of function calls used in the program
%
% Uses the following functions/m-files:
% gls.m and its subprograms
% 

function [x0,f0,l,L,istar,ncall] = initlist(fcn,data,u,v)
ncall = 0;
nloc = 5;
small = 0.1;
smaxls = 25;  
n = length(u);
x = min(max(u,0),v);	% absolutely smallest point
f = feval(fcn,data,x);
ncall = ncall + 1;
for i = 1:n
  alist = 0;
  flist = f;
  p = zeros(n,1);
  p(i) = 1;
  [alist,flist,nfls] = gls(fcn,data,u,v,x,p,alist,flist,nloc,small,smaxls);
  ncall = ncall + nfls;
  [alist1,flist1] = lspost(alist,flist);
  if isempty(find(alist1==0))
    alist1 = [alist1 0];
    flist1 = [flist1 f];
  end
  if length(alist1) < 3
    if isempty(find(alist1==alist(length(alist))))
      alist1 = [alist1 alist(length(alist))];
      flist1 = [flist1 flist(length(alist))];
    end
    if length(alist1) < 3
      if isempty(find(alist1==alist(1)))
        alist1 = [alist1 alist(1)];
        flist1 = [flist1 flist(length(alist))];
      end
      if length(alist1) < 3
        k = round((1+length(alist))/2);
        alist1 = [alist1 alist(k)];
        flist1 = [flist1 flist(k)];
      end
    end
  end
  [alist,ind] = sort(alist1);
  flist = flist1(ind);
  l(i) = find(alist == 0);
  [f1,istar(i)] = min(flist);
  L(i) = length(alist);
  x0(i,1:L(i)) = alist + x(i);
  f0(1:L(i),i) = flist';
  x(i) = x0(i,istar(i));
  f = feval(fcn,data,x);
end

