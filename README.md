# Dimension and model reduction approaches for linear Bayesian inverse problems with rank-deficient prior covariances

This MATLAB repository contains code for the numerical results of the following paper:

1. König, J., Qian E., Freitag, M. A. "[Dimension and model reduction approaches for linear Bayesian inverse problems with rank-deficient prior covariances](http://arxiv.org/abs/2506.23892)."

## Summary
The work in [1] proposes balanced truncation (BT) of a new prior-driven LTI system for model reduction in linear Bayesian inference, particularly suitable for a rank-deficient prior covariance.
Numerical examples compare the performance with BT-based Bayesian model reduction from [2] and the optimal dimension reduction from [3], both rewritten for a rank-deficient prior.
This work uses code from [2] (to be found at [https://github.com/elizqian/balancing-bayesian-inference](https://github.com/elizqian/balancing-bayesian-inference)).

## Examples
To run this code, you need the MATLAB Control System Toolbox.

To generate the plots from the paper for an incompatible prior (on the left side in Figures 1-3), run the *ISS_LR_incompat.m* script.

To generate the plots from the paper for a compatible prior (on the right side in Figures 1-3), run the *ISS_LR_compat.m* script.

## References
2. Qian, E., Tabeart, J. M., Beattie, C., Gugercin, S., Jiang, J., Kramer, P. R., and Narayan, A.
"[Model reduction for linear dynamical systems via balancing for Bayesian inference](https://link.springer.com/article/10.1007/s10915-022-01798-8)." Journal of Scientific Computing 91.29 (2022).
3. Spantini, A., Solonen, A., Cui, T., Martin, J., Tenorio, L., and Marzouk, Y. "[Optimal low-rank approximations of Bayesian linear inverse problems](https://epubs.siam.org/doi/pdf/10.1137/140977308?casa_token=CaYk5XimLkoAAAAA:-WjPu7U7kT8q3WZU66efl5X6GPylJOcnJM7XuOyy-I00LLa0vo9478Tv4BeNFoO67EwOsvl78Q)." SIAM Journal on Scientific Computing 37. 6 (2015): A2451-A2487.

### Contact
Please feel free to contact [Josie König](https://www.math.uni-potsdam.de/professuren/datenassimilation/personen/josie-koenig) with any questions about this repository or the associated paper.
