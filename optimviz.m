function funlog = optimviz(fun,x0,noise,MaxFunEvals)
%OPTIMIVIZ Optimizer visualization.
%  OPTIMVIZ runs an optimization demo with several standard optimizers on
%  the default function (Rosenbrock's function).
%  The default optimizers are define in the OPTIMDEF script.
%  The demo is also saved as an (uncompressed) animated gif. We recommend
%  to compress the animated gif before sharing.
%
%  OPTIMVIZ(FUN,X0) runs the optimization demo on the specified function
%  FUN from starting point X0 (X0 can be empty). FUN needs to have a
%  specific format (see ROSENBROCK).
%
%  OPTIMVIZ(FUN,X0,NOISE) adds i.i.d. random Gaussian noise with standard 
%  deviation NOISE to FUN (default is 0, noiseless).
%
%  OPTIMVIZ(FUN,X0,NOISE,MAXFUNEVALS) runs up to MAXFUNEVALS function 
%  evaluations (default is 100).
%
%  FUNLOG = OPTIMVIZ(...) returns the logged function for all optimizers.
%
%  See also OPTIMDEF, ROSENBROCK.

%  Author: Luigi Acerbi
%  Email: luigi.acerbi@gmail.com
%  Date: June 6, 2017
%  URL: https://github.com/lacerbi/optimviz


if nargin < 1 || isempty(fun); fun = @rosenbrock; end
if nargin < 2; x0 = []; end
if nargin < 3 || isempty(noise); noise = 0; end
if nargin < 4 || isempty(MaxFunEvals); MaxFunEvals = 100; end

seed = 0;
%ScreenPos = [1 41, 1920 958];
ScreenPos = get(0,'Screensize') + [0 40, 0 -122];
axesfontsize = 20;
fontsize = 24;

% Add subfolders to path
thispath = fileparts(mfilename('fullpath'));
addpath(genpath([thispath filesep 'utils']));

%--------------------------------------------------------------------------
% Initialize plotting and data struct

table.fun = fun;
[table.LB,table.UB,table.zlim,table.view,table.xmin,table.funname,x0_def] = table.fun();
LB = table.LB;
UB = table.UB;
if isempty(x0); x0 = x0_def; end
table.xVec = linspace(table.LB(1),table.UB(1),200);
table.yVec = linspace(table.LB(2),table.UB(2),200);
X = combvec(table.xVec,table.yVec)';
table.zMat = reshape(table.fun(X),[numel(table.xVec),numel(table.yVec)])';
table.xx = [];
table.yy = [];
table.ymin = [];
table.ymin_true = table.fun(table.xmin);
table.plotPopulation = 0;
table.funccount = 0;
table.MaxFunEvals = MaxFunEvals;
table.noise = noise;

% Added noise?
if table.noise > 0
    fun = @(x) fun(x) + table.noise*randn(size(x,1),1);
    noisy_flag = true;
else
    noisy_flag = false;
end

funlog = [];

%--------------------------------------------------------------------------
% Run optimizers

optimdef;   % Call script that defines optimizers

for iOpt = 1:numel(optims)
    rng('default');
    rng(seed);
    funlogger([],fun);
    % Run optimizer
    xs = x0;
    while 1
        optims{iOpt}{2}(@funlogger,xs,LB,UB);
        % Store results
        funlog{iOpt} = funlogger();
        if funlog{iOpt}.FuncCount >= MaxFunEvals; break; end
        
        % If execution interrupted before depleting budget, restart from
        % random point
        xs = rand(1,numel(x0)).*(UB-LB) + LB;
    end
end

%--------------------------------------------------------------------------
% Prepare figure

figure; axis off;
set(gcf,'Position',ScreenPos);
if noisy_flag; noisestr = 'noisy '; else noisestr = []; end
titlestring = ['OptimViz (' noisestr table.funname ')'];

box2 = [0.5 0.84, 0.35 0.125];  % Inner panel

switch numel(optims)
    case 1; grid = [1 1];
    case 2; grid = [1 2];
    case 3; grid = [1 3];
    case 4; grid = [2 2];
    case {5,6}; grid = [2 3];
    case {7,8}; grid = [2 4]; box2 = [0.65 0.84, 0.275 0.125];  % Inner panel
    otherwise
        error('Too many optimizers.');
end

hg = plotify(grid(1), grid(2), 'Margins', [.05 .02 .1 .15], 'Title', titlestring, 'FontSize', fontsize);
axes(hg(end)); axis off;

for iPanel = 1:numel(hg)-1
    axes(hg(iPanel));
    pos = get(hg(iPanel), 'OuterPosition');        
    pos2 = [pos(1) + box2(1)*pos(3), pos(2) + box2(2)*pos(4), box2(3)*pos(3), box2(4)*pos(4)];
    hg2(iPanel) = axes('Position', pos2);
end

axes(hg(end));
text(0.75,0.08,'https://github.com/lacerbi/bads','FontSize',axesfontsize,'HorizontalAlignment','left');
text(0.75,0.035,'https://github.com/lacerbi/optimviz','FontSize',axesfontsize,'HorizontalAlignment','left');
text(0.05,0.0525,'#useBADS','FontSize',axesfontsize,'HorizontalAlignment','left');

%--------------------------------------------------------------------------
% Plot iterations and create animated gif

splitname = strsplit(lower(table.funname));
filename = ['optimviz-' splitname{1}];
if noisy_flag; filename = [filename '-noisy']; end

for n = 1:MaxFunEvals
    for iOpt = 1:numel(funlog)    
        table.h1 = hg(iOpt);
        table.h2 = hg2(iOpt);
        table.name = optims{iOpt}{1};
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
imwrite(im,map,[filename '.gif'],'DelayTime',0,'LoopCount',inf);

end