%OPTIMDEF Define optimizers for OPTIMVIZ. 
%  Optimizers are defined by a pair of cells, the first one being the 
%  displayed name (a string), and the second a function handle that calls 
%  the optimizer, which takes as input the target function, the initial 
%  point X0, the lower bound LB and upper bound UB.

optims = [];
optims{end+1} = {'BADS', @(fun,x0,LB,UB) bads(fun,x0,LB,UB,LB,UB)};
optims{end+1} = {'fminsearch', @(fun,x0,LB,UB) fminsearch(fun,x0)};
optims{end+1} = {'fmincon', @(fun,x0,LB,UB) fmincon(fun,x0,[],[],[],[],LB,UB)};
%optims{end+1} = {'fminunc', @(fun,x0,LB,UB) fminunc(fun,x0)};
%optims{end+1} = {'patternsearch', @(fun,x0,LB,UB) patternsearch(fun,x0,[],[],[],[],LB,UB)};
optims{end+1} = {'ga', @(fun,x0,LB,UB) ga(fun,numel(x0),[],[],[],[],LB,UB,[],struct('Generations',12,'PopulationSize',25))};
optims{end+1} = {'MCS', @(fun,x0,LB,UB) mcs(func2str(fun),[],LB',UB')};
optims{end+1} = {'CMA-ES', @(fun,x0,LB,UB) cmaes_wrapper(fun,x0,LB,UB,LB,UB,struct('Seed',seed))};