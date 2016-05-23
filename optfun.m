function [y,g,flag] = optfun(x,varargin)

if nargin == 2 && isempty(x); x = varargin{1}; end

persistent table;
fontsize = 16;

if ischar(x); table = []; end

% Initialize persistent table
if isempty(table)
    table.fun = varargin{1};
    [table.LB,table.UB,table.zlim,table.view,table.xmin] = table.fun();
    table.xVec = linspace(table.LB(1),table.UB(1),200);
    table.yVec = linspace(table.LB(2),table.UB(2),200);
    X = combvec(table.xVec,table.yVec)';
    table.zMat = reshape(table.fun(X),[numel(table.xVec),numel(table.yVec)])';
    table.xx = [];
    table.yy = [];
    table.ymin = [];
    table.plotPopulation = 0;
    table.funccount = 0;
        
    if nargin > 2 && ~isempty(varargin{2}); table.noise = varargin{2}; else table.noise = 0; end
    
    close all;
    figure(1);
    set(gcf,'Position',[1 1, 1920, 1080]);    
    table.h1 = gca;
    table.h2 = axes('Position', [0.05 0.825 0.2 0.15]);
end

% Re-initialization
if ischar(x)
    table.name = x;
    y = []; g = [];
    return;
end

axes(table.h1);

% Plot function
if isstruct(x)
    if isfield(x,'swarm')   % Particle swarm plotting
        optimValues = x;
        state = varargin{1};
        
        X = optimValues.swarm;
        Y = optimValues.swarmfvals;
        
        y = 0;  % Output        
    else                    % Genetic algorithm plotting
        options = x;
        state = varargin{1};
        if nargin > 2; flag = varargin{2}; end

        X = state.Population;
        Y = state.Score;
        
        y = state;  % Output
        g = options;
        flag = 0;        
    end    
    table.plotPopulation = 1;
    table.xx = X;
    table.yy = Y;
else
    if size(x,2) == 1; x = x'; end

    [y,g] = table.fun(x);
    
    table.funccount = table.funccount + 1;
    
    % Printout for LaTeX table
    % fprintf('\\texttt{%d} & \\texttt{%.3f} & \\texttt{%.3f} & \\texttt{%.3f} \\\\\n',table.funccount,x(1),x(2),y);
    
    if isempty(table.ymin)
        table.ymin = y;
    else
        ymin = min([y,table.ymin(end)]);
        table.ymin(end+1) = ymin;
    end
    
    X = table.xx;
    Y = table.yy;

    if ~table.plotPopulation
        table.xx = [table.xx; x];
        table.yy = [table.yy; y];
    end
    
    if isa(table.noise,'function_handle')
        y = y + randn(size(y)).*table.noise(y);
    elseif table.noise > 0
        y = y + randn(size(y))*table.noise;        
    end
end

offset = 0.1;

hold off;
if isa(table.noise,'function_handle')
    Z = table.zMat + table.noise(table.zMat).*randn(size(table.zMat));
elseif table.noise > 0
    Z = table.zMat + table.noise.*randn(size(table.zMat));
else
    Z = table.zMat;
end
surf(table.xVec,table.yVec,Z,'LineStyle','none'); hold on;
scatter3(table.xmin(1),table.xmin(2),0+offset,'rd','MarkerFaceColor','r','MarkerEdgeColor','none');    % Minimum
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title(table.name,'FontSize',fontsize);
set(gca,'Xtick',table.xVec(1):table.xVec(end),'Ytick',table.yVec(1):table.yVec(end));
set(gcf,'Color','w');

if ~table.plotPopulation
    if size(table.xx,1) > 1
        plot3(table.xx(end-1:end,1),table.xx(end-1:end,2),table.yy(end-1:end,1)+offset,'-y','LineWidth',2);
    end
end
if ~isempty(X)
    scatter3(X(:,1),X(:,2),Y(:)+offset,'o','MarkerFaceColor','y','MarkerEdgeColor','none');
end
xlim([table.xVec(1),table.xVec(end)]); 
ylim([table.yVec(1),table.yVec(end)]);
zlim(table.zlim);
view(table.view);

h = text(table.UB(1),1,table.zlim(2)*5/3,['f_{evals} = ' num2str(table.funccount)],'Interpreter','TeX','FontSize',fontsize);
h = text(table.UB(1),1,table.zlim(2)*4/3,['f_{min}(x,y) = ' num2str(table.ymin(end),'%.3f')],'Interpreter','TeX','FontSize',fontsize);

axes(table.h2);
set(gca,'TickDir','out');
hold off;
semilogy(table.ymin,'-k','LineWidth',1);
box off;
xlabel('Fcn evals');
ylabel('f_{min}','Interpreter','TeX');
set(gca,'Ytick',[1e-4,1e-2,1,100]);
xlim([0,200])
ylim([1e-5,1e3]);


drawnow;
    
    
    
end

