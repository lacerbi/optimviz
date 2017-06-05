% bps_figure

N = 2^8;
fun = @rosenbrock;
fontsize = 32;
axesfontsize = 24;

rng(0);


xx = linspace(LB(1),UB(1),N);
yy = linspace(LB(2),UB(2),N);

% Compute function landscape
zz = zeros(N,N);
for j = 1:N
    zz(j,:) = fun([xx(:),yy(j)*ones(N,1)]);
end

close all;
figure;

%% Color definitions

colopt = [243,134,48]/255;
colsearch = [105,210,231]/255;

%% First panel, poll
hgrid = plotify([1 2],'Margins',[0.05 0.05],'Gutter',[0.15 0.025]);
axes(hgrid(1)); hold on; h = [];
contour(xx,yy,sqrt(zz),6);

x = [-1,0];
v1 = [0.6,0];
v2 = [0,0.3];
h(1) = plot(x(1),x(2),'ok','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k');

h(2) = plot([x(1)-v1(1),x(1)+v1(1)],[x(2)-v1(2),x(2)+v1(2)],'-dk','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2);
plot([x(1)-v2(1),x(1)+v2(1)],[x(2)-v2(2),x(2)+v2(2)],'-dk','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',2);

h(3) = plot(xmin(1),xmin(2),'pr','MarkerSize',15,'MarkerEdgeColor',colopt,'MarkerFaceColor',colopt);

xlabel('$x_1$','Interpreter','LaTeX','FontSize',fontsize);
ylabel('$x_2$','Interpreter','LaTeX','FontSize',fontsize);
set(gca,'TickDir','out','FontSize',axesfontsize);
% set(gca,'XTick',-3:3,'YTick',-3:3);
set(gca,'XTick',[],'YTick',[]);

hl = legend(h,'Current point','Polled points','Global optimum');
set(hl,'FontSize',axesfontsize);
axis square;

text(0,1.05,['Iteration 5'],'Units','Normalized','FontSize',axesfontsize);
text(1,1.05,['BADS: POLL phase'],'Units','Normalized','FontSize',fontsize,'HorizontalAlignment','right');
text(1.1,-0.05,['contours: true function'],'Units','Normalized','FontSize',axesfontsize,'HorizontalAlignment','right');

drawnow;

%% Second panel, poll
axes(hgrid(2)); hold on;

options = [];
options.FunValues.X = [x; x + v1; x - v1; x + v2; x - v2];
options.FunValues.Y = fun(options.FunValues.X);

x0 = x;
options.SkipPoll = 'on';
options.Ninit = 0;
options.SearchNtry = 20;
options.MaxFunEvals = 18;
[LB,UB,~,~,xmin] = fun();
[x,fval,exitflag,output,funValues,gpstruct] = bps(@rosenbrock,x0,LB,UB,[],[],options);

% Compute (negative) Expected Improvement (EI)
xi = combvec(xx,yy)';
[negei,~,ymu,ys,fmu,fs] = acqNegEI(xi,min(gpstruct.y(1:end-1)),gpstruct,[]);


[nn,idx] = min(negei);

negei_mat = reshape(negei,[N N])';
fi = reshape(sqrt(fmu-min(fmu(:))),[N N])';
si = reshape(fs,[N N])';
contour(xx,yy,fi,6);
hold on;
h(2) = plot(funValues.X(end,1),funValues.X(end,2),'db','MarkerSize',10,'MarkerEdgeColor',colsearch,'MarkerFaceColor',colsearch);
%h(2) = plot(xi(idx,1),xi(idx,2),'pb','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
h(1) = plot(funValues.X(1:end-1,1),funValues.X(1:end-1,2),'o','MarkerSize',6,'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',0.5*[1 1 1]);
h(3) = plot(xmin(1),xmin(2),'pr','MarkerSize',15,'MarkerEdgeColor',colopt,'MarkerFaceColor',colopt);

xlabel('$x_1$','Interpreter','LaTeX','FontSize',fontsize);
ylabel('$x_2$','Interpreter','LaTeX','FontSize',fontsize);
set(gca,'TickDir','out','FontSize',axesfontsize);
% set(gca,'XTick',-3:3,'YTick',-3:3);
set(gca,'XTick',[],'YTick',[]);

axis square;
text(0,1.05,['Iteration ' num2str(size(funValues.X,1)-1)],'Units','Normalized','FontSize',axesfontsize);
text(1,1.05,['BADS: SEARCH phase'],'Units','Normalized','FontSize',fontsize,'HorizontalAlignment','right');
text(1.1,-0.05,['contours: estimated function'],'Units','Normalized','FontSize',axesfontsize,'HorizontalAlignment','right');

% subplot(1,3,3); hold on;

xc = funValues.X(end-1,:);
rr = sqrt(bsxfun(@plus, (xx' - xc(1)).^2, (yy - xc(2)).^2)) < 1;

[~,h(4)] = contour(xx,yy,negei_mat .* rr,[1 -0.25]); hold on;
% contour(xx,yy,si,30);

%h(2) = plot(funValues.X(end,1),funValues.X(end,2),'pb','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
%h(2) = plot(xi(idx,1),xi(idx,2),'pb','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
%h(1) = plot(funValues.X(1:end-1,1),funValues.X(1:end-1,2),'o','MarkerSize',10,'MarkerEdgeColor',0.5*[1 1 1],'MarkerFaceColor',0.5*[1 1 1]);
%h(3) = plot(xmin(1),xmin(2),'pr','MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');

hl = legend(h,'Evaluated points','Next search point','Global optimum','High Expected Improvement');
set(hl,'FontSize',axesfontsize);


set(gcf,'Color','w');