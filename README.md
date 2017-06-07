## OptimViz - Optimizer visualization demo for MATLAB

This demo visualizes several MATLAB derivative-free optimizers at work on standard test functions.

The optimizers are:

- BADS (Bayesian adaptive direct search), a novel algorithm that combines a direct search approach with local Bayesian optimization ([link](https://github.com/lacerbi/bads));
- `fminsearch` (Nelder-Mead), the standard simplex method for nonlinear optimization;
- `fmincon`, a powerful first-order method;
- `ga` (genetic algorithms), a heuristic population-based method for global optimization;
- MCS (Multi-level coordinate search), an advanced method for global optimization ([link](http://www.mat.univie.ac.at/~neum/software/mcs/));
- CMA-ES (Covariance matrix adaptation - evolution strategies), a state-of-the-art method for nonconvex optimization ([link](https://www.lri.fr/~hansen/cmaesintro.html)).

We see here an example on the Rosenbrock banana function:

![demo_opt](http://luigiacerbi.com/wp-content/uploads/2016/05/demo_opt.png)

We see how the algorithms react to noise, by adding unit Gaussian noise at each function evaluation:

We see here another noiseless example on the Ackley function:


### Comments

- BADS works well on these examples, which were chosen to show how different algorithms explore the space. More generally, BADS is best for functions with a noisy or jagged landscape, and with non-negligible computational cost (see here).
- `fminsearch` is a generic optimizer which can deal with simple functions, but it should never be the main choice as there are always better alternatives.
- `fmincon` is generally superior to most optimizers (and in partcular, to `fminsearch`) on smooth functions. However, `fmincon` deals *very badly* with jagged or noisy landscapes.
- We are not aware of scenarios in which `ga` is a good off-the-shelf choice for continuous-valued optimization. It is often just barely better than random search.
- MCS can be a great optimizer, but it is somewhat idiosyncratic (it might converge very quickly to a solution).
- CMA-ES, despite the poor performance shown here, is a good optimizer *if* allowed a very large number of function evaluations.

For more details about a benchmark comparing different MATLAB optimizers on artificial and real applied problems (fitting of computational models), see the following reference:




These animated gifs can be generated via the `optimviz.m` function. 

To run some of these algorithms you will need MATLAB's [Optimization Toolbox](http://www.mathworks.com/products/optimization/) and [Global Optimization Toolbox](http://www.mathworks.com/products/global-optimization/).

