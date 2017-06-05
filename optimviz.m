function funlog = optimviz(fun,x0,noise)

if nargin < 1 || isempty(fun); fun = @rosenbrock; end
if nargin < 2 || isempty(x0); x0 = [-0.5,2.5]; end
if nargin < 3 || isempty(noise); noise = 0; end

seed = 0;
MaxFunEvals = 100;
fontsize = 24;

% Initialize function for plotting
table.fun = fun;
[table.LB,table.UB,table.zlim,table.view,table.xmin,table.funname] = table.fun();
LB = table.LB;
UB = table.UB;
table.xVec = linspace(table.LB(1),table.UB(1),200);
table.yVec = linspace(table.LB(2),table.UB(2),200);
X = combvec(table.xVec,table.yVec)';
table.zMat = reshape(table.fun(X),[numel(table.xVec),numel(table.yVec)])';
table.xx = [];
table.yy = [];
table.ymin = [];
table.plotPopulation = 0;
table.funccount = 0;
table.noise = noise;

if table.noise > 0
    fun = @(x) fun(x) + table.noise*randn(size(x,1),1);
    fun_cmaes = @(x) fun(x') + table.noise*randn(1,size(x,2));
else
    fun_cmaes = fun;
end

funlog = [];
optname = {'BADS','fminsearch','fmincon','ga','MCS','CMA-ES'};

% BADS
rng(seed);
funlogger([], fun);
bads(@funlogger,x0,LB,UB,LB,UB);
funlog{end+1} = funlogger();

% FMINSEARCH (Nelder-Mead)
rng(seed);
funlogger([], fun);
fminsearch(@funlogger,x0);
funlog{end+1} = funlogger();

% FMINCON
rng(seed);
funlogger([], fun);
fmincon(@funlogger,x0,[],[],[],[],LB,UB);
funlog{end+1} = funlogger();

% GENETIC ALGORITHM
rng(seed);
funlogger([], fun);
ga(@funlogger,2,[],[],[],[],LB,UB,[],struct('Generations',12,'PopulationSize',25));
funlog{end+1} = funlogger();

% MCS
rng(seed);
funlogger([], fun);
mcs('funlogger',[],LB',UB');
funlog{end+1} = funlogger();

% CMAES
rng(seed);
funlogger([], fun_cmaes);
options.Seed = seed; OPTS.CMA.active = 1;
cmaes_wrapper(@funlogger,x0,LB,UB,LB,UB,options); rng('default');
funlog{end+1} = funlogger();


figure; axis off;
set(gcf,'Position',[1 41, 1920, 958]);
titlestring = ['OptimViz (' table.funname ')'];

hg = plotify(2, 3, 'Margins', [.05 .02 .05 .15], 'Title', titlestring, 'FontSize', fontsize);
axes(hg(end)); axis off;

for iPanel = 1:numel(hg)-1
    axes(hg(iPanel));
    pos = get(hg(iPanel), 'OuterPosition');
    box2 = [0.5 0.84, 0.35 0.125];    
    pos2 = [pos(1) + box2(1)*pos(3), pos(2) + box2(2)*pos(4), box2(3)*pos(3), box2(4)*pos(4)];
    hg2(iPanel) = axes('Position', pos2);
end

f = getframe(gcf);

axes(hg(end));
text(0.175,0.925,'https://github.com/lacerbi/bads','FontSize',fontsize,'HorizontalAlignment','center');
text(0.825,0.925,'https://github.com/lacerbi/demo-opt','FontSize',fontsize,'HorizontalAlignment','center');

% Plot iterations
for n = 1:MaxFunEvals
    for iOpt = 1:numel(funlog)    
        table.h1 = hg(iOpt);
        table.h2 = hg2(iOpt);
        table.name = optname{iOpt};
        optimplot(funlog{iOpt},table,n);
    end
    
    f = getframe(gcf);
    if n == 1
        [im,map] = rgb2ind(f.cdata,256,'nodither');
        im(1,1,1,MaxFunEvals) = 0;
    else    
        im(:,:,1,n) = rgb2ind(f.cdata,map,'nodither');
    end
end
imwrite(im,map,'anime.gif','DelayTime',0,'LoopCount',inf)

end

%--------------------------------------------------------------------------
function y = combvec(varargin)
%COMBVEC Create all combinations of vectors.
  y = varargin{1};
  for i=2:length(varargin)
    z = varargin{i};
    y = [copy_blocked(y,size(z,2)); copy_interleaved(z,size(y,2))];
  end
end

function b = copy_blocked(m,n)
    [mr,mc] = size(m);
    b = zeros(mr,mc*n);
    ind = 1:mc;
    for i=[0:(n-1)]*mc; b(:,ind+i) = m; end
end

function b = copy_interleaved(m,n)
    [mr,mc] = size(m);
    b = zeros(mr*n,mc);
    ind = 1:mr;
    for i=[0:(n-1)]*mr; b(ind+i,:) = m; end
    b = reshape(b,mr,n*mc);
end