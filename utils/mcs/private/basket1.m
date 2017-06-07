%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% basket1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checks whether a candidate for the 'shopping basket' lies in the 
% 'domain of attraction' of a point already in the 'shopping basket'
% function [xbest,fbest,xmin,fmi,loc,flag,ncall] = basket1(fcn,data,x,f
% xmin,fmi,xbest,fbest,stop,nbasket)
% Input:
% fcn = 'fun' 	name of function fun(data,x), x an n-vector
% data		data vector (or other data structure)
% x(1:n)	candidate for the shopping basket
% f		its function value
% xmin(1:n,:)  	columns are the base vertices of the boxes in the 
%              	'shopping basket'
% fmi          	fmi(j) is the function value at xmin(:,j)
% xbest       	current best vertex
% fbest    	current best function value
% stop          stop(1) in ]0,1[:  relative error with which the known 
%		 global minimum of a test function should be found
%		 stop(2) = fglob known global minimum of a test function
%		 stop(3) = safeguard parameter for absolutely small 
%		 fglob
%		stop(1) >= 1: the program stops if the best function
%		 value has not been improved for stop(1) sweeps
%		stop(1) = 0: the user can specify a function value that
%		 should be reached
%                stop(2) = function value that is to be achieved
% nbasket	current number of points in the shopping basket
% Output:
% xbest       	current best vertex
% fbest    	current best function value
% xmin(1:n,:)	updated version of the points in the shopping basket
% fmi		their function values
% loc           = 0  candidate lies in the 'domain of attraction' of a
%		     point in the 'shopping basket 
%		= 1  otherwise
% flag		= 0  the global minimum of a test function has been 
%		     found with the required accuracy prcerr
% 		= 1  otherwise
% ncall		number of function calls used in the program

function [xbest,fbest,xmin,fmi,loc,flag,ncall] = basket1(fcn,data,x,f,xmin,fmi,xbest,fbest,stop,nbasket)
global nsweep nsweepbest 
loc = 1;
flag = 1;
ncall = 0;
if ~nbasket, return, end
for k = 1:nbasket
  dist(k) = norm(x - xmin(:,k));
end
[dist1,ind] = sort(dist);
for k = 1:nbasket
  i = ind(k);
  p = xmin(:,i) - x;
  y1 = x + 1/3*p;
  f1 = feval(fcn,data,y1);
  ncall = ncall + 1;
  if f1 <= max(fmi(i),f)
    y2 = x + 2/3*p;
    f2 = feval(fcn,data,y2);
    ncall = ncall + 1;
    if f2 <= max(f1,fmi(i))
      if f < min(min(f1,f2),fmi(i))
        fmi(i) = f;
        xmin(:,i) = x;
        if fmi(i) < fbest
          fbest = fmi(i);
          xbest = xmin(:,i);
          nsweepbest = nsweep;
          if stop(1) > 0 & stop(1) < 1
            flag = chrelerr(fbest,stop);
          elseif stop(1) == 0
            flag = chvtr(fbest,stop(2));
          end
          if ~flag,return,end
        end   
        loc = 0;break
      elseif f1 < min(min(f,f2),fmi(i))
        fmi(i) = f1;
        xmin(:,i) = y1;
        if fmi(i) < fbest
          fbest = fmi(i);
          xbest = xmin(:,i);
          nsweepbest = nsweep;
          if stop(1) > 0 & stop(1) < 1
            flag = chrelerr(fbest,stop);
          elseif stop(1) == 0
            flag = chvtr(fbest,stop(2));
          end
          if ~flag,return,end
        end   
        loc = 0;break
      elseif f2 < min(min(f,f1),fmi(i))
        fmi(i) = f2;
        xmin(:,i) = y2;
        if fmi(i) < fbest
          fbest = fmi(i);
          xbest = xmin(:,i);
          nsweepbest = nsweep;
          if stop(1) > 0 & stop(1) < 1 
            flag = chrelerr(fbest,stop);
          elseif stop(1) == 0
            flag = chvtr(fbest,stop(2));
          end
          if ~flag,return,end
        end   
        loc = 0;break
      else
        loc = 0;break
      end
    end
  end
end
