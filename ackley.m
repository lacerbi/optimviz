function varargout = ackley( x )

if nargin == 0
    varargout{1} = [-5,-5];     % LB
    varargout{2} = [5,5];       % UB
    varargout{3} = [0,20];    % Zlim
    varargout{4} = [-26.5,60];  % view
    varargout{5} = [1,1];  % xmin    
    return;
end
    varargout{1} = -20 * exp(-0.2 * sqrt( sum(x .* x, 2)/ 2)) ...
        - exp(sum(cos(2*pi*x),2) / 2) ...
        + 20 + 2.7182818284590452353602874713526625;
    varargout{2} = [];
end
