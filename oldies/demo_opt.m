%DEMO_OPT Optimization demo.
% Luigi Acerbi, May 2016

%% rosenbrock
optfun('Rosenbrock''s banana function',@rosenbrock);
optfun([1,1]);
pause 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINISTIC OBJECTIVE FUNCTION

% fminsearch
x0 = [-0.5,2.5];
optfun('fminsearch',@rosenbrock);
fminsearch(@optfun,x0);

%% fminunc
x0 = [-0.5,2.5];
optfun('fminunc',@rosenbrock);
fminunc(@optfun,x0);

%% patternsearch
x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3];
optfun('patternsearch',@rosenbrock);
patternsearch(@optfun,x0,[],[],[],[],[],[],[],struct('MaxFunEvals',300));

%% mcs
LB = [-2,-1]; UB = [2,3];
optfun('Multi-level Coordinate Search (MCS)',@rosenbrock);
try mcs('optfun',[],LB',UB'); catch display('MCS not installed.'); end

%% simulated annealing
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3];
optfun('Simulated annealing',@rosenbrock);
simulannealbnd(@optfun,x0,LB,UB,struct('MaxFunEvals',250));

%% genetic algorithm
rng('default'); rng(0); LB = [-2,-1]; UB = [2,3];
optfun('Genetic algorithm (GA)',@rosenbrock);
ga(@optfun,2,[],[],[],[],LB,UB,[],struct('Generations',12,'OutputFcns',@optfun,'PopulationSize',25));

%% particle swarm
rng('default'); rng(0); LB = [-2,-1]; UB = [2,3];
optfun('Particle swarm',@rosenbrock);
particleswarm(@optfun,2,LB,UB,optimoptions('particleswarm','OutputFcns',@(x,y) optfun(x,y),'MaxIter',20)); %struct('MaxFunEvals',250));

%% CMA-ES
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3];
optfun('CMA-ES',@rosenbrock);
cmaes_wrapper(@optfun,x0,LB,UB,LB,UB,struct('Seed',0)); rng('default');

%% Bayesian pattern search
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3]; warning off;
%optfun('Bayesian pattern search (BPS)',@rosenbrock);
% try bps(@optfun,x0,LB,UB); catch display('BPS not installed.'); end
optfun('Bayesian adaptive direct search (BADS)',@rosenbrock);
try bads(@optfun,x0,LB,UB,LB,UB,[],struct('TolFun',1e-3)); catch display('BADS not installed.'); end

%% IMGPO
% rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3];
% optfun('IMGPO',@rosenbrock);
% try IMGPO_default_run(@optfun,numel(LB),[LB;UB]',300,0,0,0,0); catch display('IMGPO not installed.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOISY OBJECTIVE FUNCTION

%% noisy fminsearch
rng('default'); rng(0); x0 = [-0.5,2.5]; noise = @(y) 1 + 0.5*log(1+y);
optfun('fminsearch (noisy)',@rosenbrock,noise);
fminsearch(@optfun,x0,struct('MaxFunEvals',150));

%% noisy fminunc
rng('default'); rng(0); x0 = [-0.5,2.5]; noise = @(y) 1 + 0.5*log(1+y);
optfun('fminunc (noisy)',@rosenbrock,noise);
fminunc(@optfun,x0);

%% noisy patternsearch
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3]; noise = @(y) 1 + 0.5*log(1+y);
optfun('patternsearch (noisy)',@rosenbrock,noise);
patternsearch(@optfun,x0,[],[],[],[],[],[],[],struct('MaxFunEvals',300));

%% noisy mcs
rng('default'); rng(0); LB = [-2,-1]; UB = [2,3]; noise =@(y) 1 + 0.5*log(1+y);
optfun('MCS (noisy)',@rosenbrock,noise);
mcs('optfun',[],LB',UB');

%% noisy genetic algorithm
rng('default'); rng(0); LB = [-2,-1]; UB = [2,3]; noise = @(y) 1 + 0.5*log(1+y);
optfun('GA (noisy)',@rosenbrock,noise);
ga(@optfun,2,[],[],[],[],LB,UB,[],struct('Generations',12,'OutputFcns',@optfun,'PopulationSize',25));

%% noisy CMA-ES
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-4,-4]; UB = [5,5]; noise = @(y) 1 + 0.5*log(1+y);
optfun('CMA-ES',@rosenbrock,noise);
cmaes_wrapper(@optfun,x0,LB,UB,LB,UB,struct('Seed',0,'Noise',struct('on',1))); rng('default');

%% noisy Bayesian pattern search
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-2,-1]; UB = [2,3]; noise = @(y) 1 + 0.5*log(1+y);
%optfun('BPS (noisy)',@rosenbrock,noise);
%try bps(@optfun,x0,LB,UB,[],[],struct('UncertaintyHandling','on')); catch display('BPS not installed.'); end
optfun('BADS (noisy)',@rosenbrock,noise);
try bads(@optfun,x0,LB,UB,LB,UB,[],struct('UncertaintyHandling','on')); catch display('BADS not installed.'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIMODAL OBJECTIVE FUNCTION

%% genetic algorithm
rng('default'); rng(0); LB = [-4,-4]; UB = [5,5];
optfun('Genetic algorithm (GA)',@ackley);
ga(@optfun,2,[],[],[],[],LB,UB,[],struct('Generations',12,'OutputFcns',@optfun,'PopulationSize',25));

%% CMA-ES
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-4,-4]; UB = [5,5];
optfun('CMA-ES',@ackley);
cmaes_wrapper(@optfun,x0,LB,UB,LB,UB,struct('Seed',0)); rng('default');

%% mcs
rng('default'); rng(0); LB = [-4,-4]; UB = [5,5];
optfun('Multi-level Coordinate Search (MCS)',@ackley);
try mcs('optfun',[],LB',UB'); catch display('MCS not installed.'); end

%% Bayesian pattern search
rng('default'); rng(0); x0 = [-0.5,2.5]; LB = [-4,-4]; UB = [5,5]; warning off;
optfun('Bayesian pattern search (BPS)',@ackley);
try bps(@optfun,x0,LB,UB); catch display('BPS not installed.'); end