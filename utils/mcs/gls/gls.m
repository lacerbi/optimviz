

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gls.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [alist,flist,nf]=gls(func,data,xl,xu,x,p,alist,flist,nloc,small,smax,prt)
% local and global line search for noiseless functions
% 
% function [alist,flist,nf]=gls(prt) 	
% debug mode: repeats last line search with new print level 
%
% func:		name for function f=f(data,x)
% data: 	data vector for function call
% xl,xu:	box [xl,xu] for x
% x: 		starting point
% p:		search direction
% alist,flist:	list of known steps and corresponding function values 
%		(at output in order of increasing alp)
% nloc:		saturate nloc best local optima (default 1)
% small:	saturation tolerance (default 0.1)
% smax:		approximate limit of list size (default 10)
% prt:		print level (default 0: nothing plotted or printed)
% 		1: plot progress along the line, and
% 		2: print some things
% 		3: print more things
% nf:		number of function evaluations used
%
% Note: This is a very slow implementation since lots of sort and
%	reshape are used for simplicity of programming
%	but the number of function values used is low
%
% Note: subroutines without arguments 
%	return with a sorted alist and an updated s=size(alist,2);
% 
function [alist,flist,nf]=gls(func,data,xl,xu,x,p,alist,flist,nloc,small,smax,prt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% tuning parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
short=0.381966;			% fraction for splitting intervals
% golden section fraction is (3-sqrt(5))/2= 0.38196601125011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% set debug info %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save input for analysis of problems in debug mode
if nargin==1, 
  disp(' ');
  disp('******** debug mode; old data are used ********')
  disp(' ');
  prt1=func;
  load glsinput;
  prt=prt1;
else 
  if nargin<9, nloc=1; end;
  if nargin<10, small=0.1; end;
  if nargin<11, smax=10; end;
  if nargin<12, prt=0; end;
  % save glsinput;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if prt>1, 
  disp('********************** gls *********************');
  format short; 
end;

% save information for nf computation and extrapolation decision
sinit=size(alist,2);		% initial list size


% get 5 starting points (needed for lslocal)
bend=0;lsrange;			% find range of useful alp
lsinit; 			% 2 points needed for lspar and lsnew
nf=s-sinit;			% number of function values used

while s<min(5,smax), 
  if nloc==1,
    lspar;			% parabolic interpolation step
    % when s==3 we haven't done a true parabolic step 
    % and may appear monotonic without being so!
    if s>3 & monotone & (abest==amin | abest==amax), 
      if prt>1, disp('return since monotone'); end;
      nf=s-sinit;		% number of function values used
      lsdraw; return;
    end;
  else 
    lsnew;			% extrapolation or split
  end;
end;

saturated=0;			% is reset in lsquart
% shape detection phase
if nmin==1, 
  if monotone & (abest==amin | abest==amax), 
    if prt>1, disp('return since monotone'); end;
    nf=s-sinit;			% number of function values used
    lsdraw; return;
  end;
  if s==5, lsquart; end;	% try quartic interpolation step
  lsdescent;  			% check descent condition 		
  lsconvex; 			% check convexity	
  if convex,
    if prt>1, disp('return since convex'); end;
    nf=s-sinit;			% number of function values used
    lsdraw; return;
  end;
end;
sold=0;
% refinement phase
while 1,
  lsdraw;
  if prt>1, disp('***** new refinement iteration *****'); end;
  lsdescent;  			% check descent condition 		
  lssat;			% check saturation
  if saturated | s==sold | s>=smax, 
    if saturated & prt>1, disp('return since saturated'); end;
    if s==sold   & prt>1, disp('return since s==sold');   end;
    if s>=smax   & prt>1, disp('return since s>=smax');   end;
    lsdraw;break; 
  end;
  sold=s; nminold=nmin; if prt>1, nmin, end;
  if ~saturated & nloc>1,
    lssep; 			% separate close minimizers
  end;
  lslocal;			% local interpolation step
  if nmin>nminold, saturated=0; end;
end;

% get output information
nf=s-sinit;			% number of function values used



