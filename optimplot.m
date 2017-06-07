function optimplot(funlog,table,n)
%OPTIMPLOT Plot optimization panel for OPTIMVIZ.
%
%  See also OPTIMVIZ.

fontsize = 24;
axesfontsize = 18;

noisy_flag = ~isempty(table.noise) && table.noise > 0;

axes(table.h1);
        
    % Printout for LaTeX table
    % fprintf('\\texttt{%d} & \\texttt{%.3f} & \\texttt{%.3f} & \\texttt{%.3f} \\\\\n',table.funccount,x(1),x(2),y);
        
X = funlog.X(1:n,:);

% Re-evaluate Y with true function (assumed vectorized)
Y = table.fun(X);
E = Y - table.ymin_true;

% Build running error
ymin = zeros(n,1);
for i = 1:n; ymin(i) = min(E(1:i)); end

%     if ~table.plotPopulation
%         table.xx = [table.xx; x];
%         table.yy = [table.yy; y];
%     end
    
%     if isa(table.noise,'function_handle')
%         y = y + randn(size(y)).*table.noise(y);
%     elseif table.noise > 0
%         y = y + randn(size(y))*table.noise;        
%     end

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
scatter3(table.xmin(1),table.xmin(2),0+offset,100,'rd','MarkerFaceColor','r','MarkerEdgeColor','none');    % Minimum
xlabel('x','FontSize',axesfontsize);
ylabel('y','FontSize',axesfontsize);
zlabel('f(x,y)','FontSize',axesfontsize);
% title(table.name,'FontSize',fontsize);
set(gca,'Xtick',table.xVec(1):table.xVec(end),'Ytick',table.yVec(1):table.yVec(end));
set(gcf,'Color','w');
set(gca,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);


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
grid off; axis off;
view(table.view);

h = text(0.15,0.965,table.name,'FontSize',fontsize,'HorizontalAlignment','center','Units','Normalized');

% 
% h = text(table.UB(1),1,table.zlim(2)*5/3,['f_{evals} = ' num2str(table.funccount)],'Interpreter','TeX','FontSize',fontsize);
% h = text(table.UB(1),1,table.zlim(2)*4/3,['f_{min}(x,y) = ' num2str(table.ymin(end),'%.3f')],'Interpreter','TeX','FontSize',fontsize);

axes(table.h2);
set(gca,'TickDir','out');
hold off;
semilogy(ymin,'-k','LineWidth',3);
box off;
xlabel('Fcn. evals.','FontSize',axesfontsize);
ylabel('Error','FontSize',axesfontsize);
xlim([0,100]);

if noisy_flag
    ylim([1e-2,1e3]);
    ytick = [1e-2,1,1e2];
    yticklabel = {'0.01','1','100'};
else
    ylim([1e-3,1e3]);
    ytick = [1e-3,1,1e3];
    yticklabel = {'0.001','1','1000'};
end
set(gca,'XTick',[0:50:300],'Ytick',ytick,'YTickLabel',yticklabel,'TickDir','out','FontSize',axesfontsize-4);

drawnow;
    
end
