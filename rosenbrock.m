function varargout = rosenbrock(x)
%ROSENBROCK Rosenbrock's 'banana' function in 2D

if nargin == 0
    varargout{1} = [-2,-1];                 % LB
    varargout{2} = [2,3];                   % UB
    varargout{3} = [0,3000];                % Zlim
    varargout{4} = [-26.5,60];              % view
    varargout{5} = [1,1];                   % xmin
    varargout{6} = 'Rosenbrock function';   % name
    varargout{7} = [-0.5,2.5];              % default x0
    return;
end

if size(x,2) == 1; x = x'; end

varargout{1} = 100*(x(:,2)-x(:,1).^2).^2+(1-x(:,1)).^2;
if nargout > 1
    varargout{2} = [-400*(x(:,2)-x(:,1).^2).*x(:,1)-2*(1-x(:,1)), 200*(x(:,2)-x(:,1).^2)];
end

end