%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% defaults.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [u,v,fglob] = defaults(fcn)
% sets the default original boxes [u,v] and the optimal function value
% fglob for the Jones test functions fcn in {'bra','cam','gpr','hm3',
% 'hm6','s10','sh5','sh7','shu'};
% in addition the global variables fglob, xglob and nglob are defined
function [u,v,fglob] = defaults(fcn)
global nglob xglob
% fglob        global minimum of fcn
% nglob        number of global minimizers in the box [u,v]
% xglob(1:n,1:nglob)  xglob(:,i), i = 1:nglob, are the global
%              minimizers of a test function fcn in [u,v]
if fcn == 'gpr'      % Goldstein-Price
  u = [-2; -2]; 
  v = [2; 2]; 
  fglob = 3;
  xglob = [0.; -1.];
  nglob = 1;
elseif fcn == 'bra'  % Branin
  u = [-5; 0];
  v = [10; 15]; 
  fglob = 0.397887357729739;
  xglob = [9.42477796   -3.14159265  3.14159265; 
           2.47499998   12.27500000  2.27500000];
  nglob = 3;  
elseif fcn == 'cam'  % Six-hump camel
  u = [-3; -2];
  v = [3; 2]; 
  fglob = -1.0316284535;
  xglob = [ 0.08984201  -0.08984201;
           -0.71265640   0.71265640];
  nglob = 2;
elseif fcn == 'shu'  % Shubert
  n = 2;
  u = [-10; -10];
  v = [10; 10];
  fglob = -186.730908831024;
  xglob = [
-7.08350658  5.48286415  4.85805691  4.85805691 -7.08350658 -7.70831382 -1.42512845 -0.80032121 -1.42512844 -7.08350639 -7.70831354  5.48286415  5.48286415  4.85805691 -7.70831354 -0.80032121 -1.42512845 -0.80032121; 
 4.85805691  4.85805681 -7.08350658  5.48286415 -7.70831382 -7.08350658 -0.80032121 -1.42512845 -7.08350639 -1.42512844  5.48286415 -7.70831354  4.85805691  5.48286415 -0.80032121 -7.70831354 -0.80032121 -1.42512845];
  nglob = 18;
elseif fcn == 'sh5'  % Shekel 5
  u = [0; 0; 0; 0];
  v = [10; 10; 10; 10]; 
  fglob = -10.1531996790582;
  xglob = [4; 4; 4; 4];
  nglob = 1;
elseif fcn == 'sh7'  % Shekel 7
  u = [0; 0; 0; 0];
  v = [10; 10; 10; 10]; 
  fglob = -10.4029405668187;
  xglob = [4; 4; 4; 4];
  nglob = 1;
elseif fcn == 's10'  % Shekel 10
  u = [0; 0; 0; 0];
  v = [10; 10; 10; 10]; 
  fglob = -10.5364098166920;
  xglob = [4; 4; 4; 4];
  nglob = 1;
elseif fcn == 'hm3'  % Hartman 3
  u = [0; 0; 0];
  v = [1; 1; 1]; 
  fglob = -3.86278214782076;
  xglob = [0.1; 0.55592003; 0.85218259];
  nglob = 1;
elseif fcn == 'hm6'  % Hartman 6
  u = [0; 0; 0; 0; 0; 0];
  v = [1; 1; 1; 1; 1; 1]; 
  fglob = -3.32236801141551;
  xglob = [0.20168952;  0.15001069;  0.47687398;  0.27533243;  0.31165162;  0.65730054];
  nglob = 1;
end
