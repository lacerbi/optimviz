

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% func.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function f=func(data,x);
% for testing gls.m, data is unused and x is only a single argument
%
function f=func(data,x);


cas=4;

if cas==1,
  f=(((3*x-1).*x-2).*x-1);
elseif cas==2,
  f=(3*x-1).^3;
elseif cas==3,
  f=(3*x-1)^2;
elseif cas==4,
  f=(x+1).*sin(10*x+1);
end;



