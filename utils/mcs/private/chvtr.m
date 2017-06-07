%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chvtr.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function chvtr(f,vtr)
% checks whether a required value to reach has already been reached; in
% that case flag is set to 0, otherwise it is set to 1
% Input:
% f	function value to be checked
% vtr	value to reach
% Output:
% flag	= 0  vtr has been reached
%	= 1  otherwise

function flag = chvtr(f,vtr)
if f <= vtr
  flag = 0;
else
  flag = 1;
end
