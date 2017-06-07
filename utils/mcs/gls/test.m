
% test program for testing gls
% adapt function in testfun.m
% 

clear;clear mex;clf;
prt=2;		% print level
xx=[-9,2];	% interval for searching minimum

func='testfun';
data=[];	% no data in test function
xl=xx(1);xu=xx(2);
x=(xl+xu)/2;
p=1;
alist=[];flist=[];
nloc=3;small=0.1;smax=10;

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

[alist,flist,nf]=gls(func,data,xl,xu,x,p,alist,flist,nloc,small,smax,prt)
