%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initbox.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ipar,level,ichild,f,isplit,p,xbest,fbest]=initbox(x0,f0,l,
% L,istar,u,v,prt)
% generates the boxes in the initialization procedure
% For each coordinate i we compute the quadratic interpolant to any
% three consecutive (x0(i,j),f0(j,i)). xl and xu are the minimum and
% maximum of the quadratic interpolant, respectively. If x0(i,istar(i)) 
% belongs to two boxes, the one containing the minimizer of the 
% quadratic model in coordinate i is taken as the current box for the
% next coordinate. All boxes generated in the initialization procedure
% get the level of their parent box increased by one.
%
% Input:
% x0(1:n,1:max(L))   coordinates used in the initialization list
% f0(1:max(L),1:n)  function values appertaining to the initialization
%             list; f0(j,i) = function value at x = (x0(1,istar(1)),...,
%             x0(i-1,istar(i-1)),x0(i,j),x0(i+1,l(i+1)),...,x0(n,l(n)))
% l(1:n)      pointer to the initial point x^0 = (x0(1,l(1)),...,
%             x0(n,l(n)) of the initialization list
% L(1:n)      the init. list contains L(i) values for coordinate i
% istar(1:n)  pointer to the final best point x^* of the init. list
%             x^* = (x0(1,istar(1)),...,x0(n,istar(n))) 
% [u,v]       initial box
% Output:
% ipar        vector containing the labels of the parents of the newly
%             created boxes
% level       vector containing their levels
% ichild      the absolute values of this vector specify which child 
%             they are; ichild(j) < 0 for the boxes generated in the 
%             initialization procedure
% f(1:2,:)    f(1,j) contains the base vertex function value of box j
% isplit      vector containing the splitting indices of the boxes that
%             have already been split
% p(1:n)      ranking of variability; the ith component has the p(i)th 
%             highest estimated variability 
% xbest       best vertex after the initialization procedure
% fbest       best function value after the initialization procedure
% prt	      print level
%
% Uses the following functions/m-files:
% genbox.m
% polint.m
% quadmin.m
% quadpol.m 
% split1.m
% updtoptl.m 

function [ipar,level,ichild,f,isplit,p,xbest,fbest]=initbox(x0,f0,l,L,istar,u,v,prt)
global nboxes nglob xglob
% nboxes      counter for boxes not in the base list
% nglob       number of global minimizers of a test function in [u,v]
% xglob(1:n,1:nglob) x(:,i), i=1:nglob, are the global minimizers of a 
% test function in [u,v]
n = length(u);
% parameter values appertaining to box 1
ipar(1) = 0;
level(1) = 1;  
ichild(1) = 1;
f(1,1) = f0(l(1),1);
par = 1;
if prt > 1,
  iopt = 1:nglob;
end
for i = 1:n 
  isplit(par) = -i;  
  % boxes split in the init. procedure get a negative splitting index
  nchild = 0;
  if x0(i,1) > u(i)
    nboxes = nboxes + 1;
    nchild = nchild + 1;
    [ipar(nboxes),level(nboxes),ichild(nboxes),f(1,nboxes)] = genbox(par,level(par)+1,-nchild,f0(1,i)); 
    if prt > 1,
      updtoptl(i,u(i),x0(i,1),iopt,level(par)+1,f0(1,i));
    end
  end
  if L(i) == 3
    v1 = v(i);
  else
    v1 = x0(i,3);
  end
  d = polint(x0(i,1:3),f0(1:3,i));
  xl = quadmin(u(i),v1,d,x0(i,1:3));
  fl = quadpol(xl,d,x0(i,1:3));
  xu = quadmin(u(i),v1,-d,x0(i,1:3));
  fu = quadpol(xu,d,x0(i,1:3));
  if istar(i) == 1
    if xl < x0(i,1)
      par1 = nboxes;  % label of the current box for the next coordinate
      j1 = 0;    % auxiliary index  
    else
      par1 = nboxes + 1;
      j1 = 2;
    end  
  end
  for j = 1:L(i)-1
    nboxes = nboxes + 1;
    nchild = nchild + 1;
    if f0(j,i) <= f0(j+1,i)
      s = 1;
    else
      s = 2;
    end
    [ipar(nboxes),level(nboxes),ichild(nboxes),f(1,nboxes)] = genbox(par,level(par)+s,-nchild,f0(j,i));
    if prt
      splval = split1(x0(i,j),x0(i,j+1),f0(j,i),f0(j+1,i));
    end
    if prt > 1,
      updtoptl(i,x0(i,j),splval,iopt,level(par)+1,f0(j,i));
    end
    if j >= 2
      if istar(i) == j
        if xl <= x0(i,j)
          par1 = nboxes - 1;
          j1 = j - 1;
        else
          par1 = nboxes;
          j1 = j + 1;
        end
      end
      if j <= L(i)-2
        d = polint(x0(i,j:j+2),f0(j:j+2,i));
        if j < L(i)-2
          u1 = x0(i,j+2);
        else
          u1 = v(i);
        end
        xl = quadmin(x0(i,j),u1,d,x0(i,j:j+2));
        fl = min(quadpol(xl,d,x0(i,j:j+2)),fl);
        xu = quadmin(x0(i,j),u1,-d,x0(i,j:j+2));
        fu = max(quadpol(xu,d,x0(i,j:j+2)),fu);
      end
    end
    nboxes = nboxes + 1;
    nchild = nchild + 1; 
    [ipar(nboxes),level(nboxes),ichild(nboxes),f(1,nboxes)] = genbox(par,level(par)+3-s,-nchild,f0(j+1,i));
    if prt > 1,
      updtoptl(i,splval,x0(i,j+1),iopt,level(par)+1,f0(j+1,i));
    end
  end
  if x0(i,L(i)) < v(i)
    nboxes = nboxes + 1;
    nchild = nchild + 1;
    [ipar(nboxes),level(nboxes),ichild(nboxes),f(1,nboxes)] = genbox(par,level(par)+1,-nchild,f0(L(i),i));
    if prt > 1,
      updtoptl(i,x0(i,L(i)),v(i),iopt,level(par)+1,f0(L(i),i));
    end
  end
  if istar(i) == L(i)
    if x0(i,L(i)) < v(i)
      if xl <= x0(i,L(i))
        par1 = nboxes - 1;
        j1 = L(i) - 1;
      else
        par1 = nboxes;
        j1 = L(i) + 1;
      end
    else
      par1 = nboxes;
      j1 = L(i) - 1;
    end
  end
  var(i) = fu - fl;  
  % the quadratic model is taken as a crude measure of the variability 
  % in the ith component
  level(par) = 0;    % box is marked as split
  par = par1;
  if j1 == 0
    splval = u(i);
  elseif j1 == L(i) + 1;
    splval = v(i);
  else
    if j1 < istar(i)
      splval = split1(x0(i,j1),x0(i,istar(i)),f0(j1,i),f0(istar(i),i));
    else
      splval = split1(x0(i,istar(i)),x0(i,j1),f0(istar(i),i),f0(j1,i));
    end
  end     
  if prt > 1 & i <= n - 1
    iopt1 = [];  
    for j = 1:length(iopt)
      if min(splval,x0(i,istar(i))) <= xglob(i,iopt(j)) & xglob(i,iopt(j)) <= max(splval,x0(i,istar(i))) 
        iopt1 = [iopt1, iopt(j)];
      end
    end  
    iopt = iopt1;      
  end
end
fbest = f0(istar(n),n);	  % best function value after the init. procedure
for i = 1:n
  [var0,p(i)] = max(var);
  var(p(i)) = -1;
  xbest(i) = x0(i,istar(i));  % best point after the init. procedure
end
