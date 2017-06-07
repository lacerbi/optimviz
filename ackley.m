function varargout = ackley(x)
%ACKLEY Ackley's function

if nargin == 0
    varargout{1} = [-8,-8];     % LB
    varargout{2} = [2,2];       % UB
    varargout{3} = [0,20];      % Zlim
    varargout{4} = [-26.5,60];  % view
    varargout{5} = [0,0];       % xmin
    varargout{6} = 'Ackley function';    % name
    varargout{7} = [-6,-5];     % default x0    
    return;
end
    varargout{1} = -20 * exp(-0.2 * sqrt( sum(x .* x, 2)/ 2)) ...
        - exp(sum(cos(2*pi*x),2) / 2) ...
        + 20 + 2.7182818284590452353602874713526625;
    varargout{2} = [];
end
