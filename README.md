# OptimViz - Optimizer visualization demo for MATLAB

This demo visualizes several MATLAB derivative-free optimizers at work on standard test functions. This is purely for demonstration purposes. For a proper benchmark of different MATLAB optimizers, see [[1](https://github.com/lacerbi/optimviz#reference)].

[Follow me on Twitter](https://twitter.com/AcerbiLuigi) for updates about other projects I am involved with, or drop me an email at  <luigi.acerbi@nyu.edu> to talk about computational modeling, optimization, and MCMC.

I have been giving seminars and tutorials on optimization, model fitting, and model comparison. If you are interested, see [my webpage](http://luigiacerbi.com/tutorials/).

## Optimizers

The optimization algorithms visualized here are:

- BADS (*Bayesian adaptive direct search*), a novel algorithm that combines a direct search approach with local Bayesian optimization ([link](https://github.com/lacerbi/bads));
- `fminsearch` (Nelder-Mead), the standard simplex method for nonlinear optimization;
- `fmincon`, a powerful method for constrained optimization based on numerical approximation of the gradient;
- `ga` (genetic algorithms), a heuristic population-based method for global optimization;
- MCS (*Multi-level coordinate search*), an advanced method for global optimization ([link](http://www.mat.univie.ac.at/~neum/software/mcs/));
- CMA-ES (*Covariance matrix adaptation - evolution strategies*), a state-of-the-art method for nonconvex optimization ([link](https://www.lri.fr/~hansen/cmaesintro.html)).

## Examples

We see here an example on the Rosenbrock banana function:

![demo_opt](https://github.com/lacerbi/optimviz/blob/master/gifs/optimviz-rosenbrock.gif)

We see how the algorithms react to noise, by adding unit Gaussian noise at each function evaluation:

![demo_opt](https://github.com/lacerbi/optimviz/blob/master/gifs/optimviz-rosenbrock-noisy.gif)

We see here another noiseless example on the Ackley function:

![demo_opt](https://github.com/lacerbi/optimviz/blob/master/gifs/optimviz-ackley.gif)


## Comments

- BADS works well on these examples, which were chosen to show how different algorithms explore the space. More generally, BADS is best for functions with a noisy or jagged landscape, and with non-negligible computational cost (see [here](https://github.com/lacerbi/bads/wiki#which-kind-of-problems-is-bads-suited-for)). BADS is available as a ready-to-use MATLAB toolbox [here](https://github.com/lacerbi/bads).
- `fminsearch` is a generic optimizer which can deal with simple functions, but it should never be the main choice as there are always better alternatives.
- `fmincon` is generally superior to most optimizers (and in partcular, to `fminsearch`) on smooth functions. However, `fmincon` deals *very badly* with jagged or noisy landscapes.
- We are not aware of scenarios in which `ga` is a good off-the-shelf choice for continuous-valued optimization. It is often just barely better than random search.
- MCS can be a great optimizer, but it is somewhat idiosyncratic (it might converge very quickly to a solution).
- CMA-ES, despite the poor performance shown here, is a good optimizer *if* allowed a very large number of function evaluations.

## Code

These animated gifs can be generated via the `optimviz.m` function. You can easily test different optimizers and add other functions.

The generated animated gifs are uncompressed. We recommend to compress them before using them in any form (e.g., via some [online tool](https://ezgif.com/optimize)).

To run some of these algorithms you will need MATLAB's [Optimization Toolbox](http://www.mathworks.com/products/optimization/) and [Global Optimization Toolbox](http://www.mathworks.com/products/global-optimization/).

## References

For more details about the benchmark comparing different MATLAB optimizers on artificial and real applied problems (fitting of computational models), see the following reference:

1. Acerbi, L. & Ma, W. J. (2017). Practical Bayesian Optimization for Model Fitting with Bayesian Adaptive Direct Search. In *Advances in Neural Information Processing Systems 30*, pages 1834-1844. ([link](https://papers.nips.cc/paper/6780-practical-bayesian-optimization-for-model-fitting-with-bayesian-adaptive-direct-search), [arXiv preprint](https://arxiv.org/abs/1705.04405))

For more info about my work in computational neuroscience and machine learning, follow me on Twitter: https://twitter.com/AcerbiLuigi


### License

*OptimViz* is released under the terms of the [GNU General Public License v3.0](https://github.com/lacerbi/optimviz/blob/master/LICENSE.txt).

