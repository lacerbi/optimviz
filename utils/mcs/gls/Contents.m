

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GLS - univariate local or global optimization           %%%%%%%%
%%%%%%% source: http://www.mat.univie.ac.at/~neum/software/ls %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% implemented as optimization in R^n along a ray x+alp*p (alp>=0)
%
% In particular, it may be used for coordinate searches by choosing 
% a unit basis vector for the search direction p
%
% For further comments see the Readme
%
% gls.m			1D line search, main routine
%			rerun with gls(prt) for detailed output
% glsinput.mat		saved last gls input parameters
% lspost.m		postprocessing of gls.m output 
% 
%
% gls.m calls the following routines:
% eqstr.m		checks equality of two strings
% lsconvex.m		check convexity
% lsdescent.m		check descent condition
% lsdraw.m		draw function along line
% lsguard.m		safeguard a step
% lsinit.m		find first two points
% lslocal.m		local refinement
% lsnew.m		extrapolation or split
% lspar.m		local parabolic interpolation step
% lsquart.m		quartic interpolation step
% lsrange.m		compute good range for step
% lssat.m		check saturation
% lssep.m		separate close local minimizers
% lssort.m		sort list by increasing step
% lssplit.m		find split of interval towards better point
% quartic.m		quartic evaluation
%
% test gls.m with
% test.m		test program
% testfun.m		test function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gls.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [alist,flist,nf]=gls(func,data,xl,xu,x,p,alist,flist,nloc,small,smax,prt)
% local and global line search for noiseless functions
% 
% function [alist,flist,nf]=gls(prt) 
% debug mode: repeats last line search with new print level 
%
% func:         name for function f=f(data,x)
% data:         data vector for function call
% xl,xu:        box [xl,xu] for x
% x:            starting point
% p:            search direction
% alist,flist:  list of known steps and corresponding function values 
%               (at output in order of increasing alp)
% nloc:         saturate nloc best local optima (default 1)
% small:        saturation tolerance (default 0.1)
% smax:         approximate limit of list size (default 10)
% prt:          print level (default 0: nothing plotted or printed)
%               1: plot progress along the line, and
%               2: print some things
%               3: print more things
% nf:           number of function evaluations used
%
% Note: This is a very slow implementation since lots of sort and
%       reshape are used for simplicity of programming
%       but the number of function values used is low
%
% Note: subroutines without arguments 
%       return with a sorted alist and an updated s=size(alist,2);
% 

