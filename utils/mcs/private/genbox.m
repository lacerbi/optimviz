%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% genbox.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ipar,level,ichild,f] = genbox(par,level0,nchild,f0)
% generates a box with parent box # ipar = par, level = level0, ichild 
% = nchild and base vertex function value f = f0 

function [ipar,level,ichild,f] = genbox(par,level0,nchild,f0)
ipar = par;
level = level0;
ichild = nchild;
f = f0;
