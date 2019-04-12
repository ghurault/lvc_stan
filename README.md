# lvc_stan
Fitting a Lotka-Volterra competition model using Stan

The data is available in two csv files (one for each set of experiment) and is from:
B. Y. G. F. Gause,
“EXPERIMENTAL STUDIES ON THE STRUGGLE FOR EXISTENCE I . MIXED POPULATION OF TWO SPECIES OF YEAST”,
J. Exp. Biol., vol. 9, no. 4, pp. 389–402, 1932.

Two models are proposed, one with multiplicative error, the other with additive error.
The two models notably implement a Cauchy regularisation on the alpha matrix (matrix of competition rates). Possibility to change the distribution to a Gaussian (cf. ridge regularisation), Laplace (cf. lasso regularisation) or to the Horseshoe distribution.

Coefficients estimates for the different models (and those estimated in the paper) are avaible in results.jpg
