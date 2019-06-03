# Lotka-Volterra competition model in Stan

## Data

The data is available in two csv files in the `Data` folder (one for each set of experiment) and was published in [Gause (1932)](http://jeb.biologists.org/content/jexbio/9/4/389.full.pdf) to study the competition of two species of yeast: saccharomyces and schixosacharomyces.

The data was collected in 3 differents experimental settings, each of them repeated twice (experiment A and B):
- Saccharomyces alone (single species setting)
- Schixosacharomyces alone (single species setting)
- Saccharomyces and Schixosacharomyces together (mixed species setting)

![data](Data/data.jpg)

## Model

The evolution of the population of yeasts is modelled with the competitive Lokta-Volterra equations.

### Lotka-Volterra competition model

For a single species p:

<img src="https://latex.codecogs.com/svg.latex?\frac{dp}{dt}&space;=&space;p&space;(r&space;-&space;\alpha&space;p)" title="\frac{dp}{dt} = p (r - \alpha p)" />

<img src="https://latex.codecogs.com/svg.latex?p(t)&space;=&space;\frac{k}{1&space;&plus;&space;C&space;e^{-&space;r&space;t}}" title="p(t) = \frac{k}{1 + C e^{- r t}}" />

- r is the growth parameter
- alpha is the (self) competition rate
- k = r / alpha is the carrying capacity (steady state)
- C is a constant determined by initial conditions

For two species:

<img src="https://latex.codecogs.com/svg.latex?\begin{cases}\frac{dp_1}{dt}&space;=&space;p_1&space;(r_1&space;-&space;\alpha_{1,&space;1}&space;p_1&space;-&space;\alpha_{1,2}&space;p_2)&space;\\\frac{dp_2}{dt}&space;=&space;p_2&space;(r_2&space;-&space;\alpha_{2,&space;1}&space;p_1&space;-&space;\alpha_{2,&space;2}&space;p_2)\end{cases}" title="\begin{cases}\frac{dp_1}{dt} = p_1 (r_1 - \alpha_{1, 1} p_1 - \alpha_{1,2} p_2) \\\frac{dp_2}{dt} = p_2 (r_2 - \alpha_{2, 1} p_1 - \alpha_{2, 2} p_2)\end{cases}" />

More generally, for N species (`*` is the element-wise multiplication):

<img src="https://latex.codecogs.com/svg.latex?\frac{d\boldsymbol{p}}{dt}&space;=&space;\boldsymbol{p}&space;*&space;(\boldsymbol{r}&space;-&space;A&space;\boldsymbol{p})" title="\frac{d\boldsymbol{p}}{dt} = \boldsymbol{p} * (\boldsymbol{r} - A \boldsymbol{p})" />

### Statistical model

The solutions of the ODEs can expressed as function of the initial conditions f<sub>0</sub>, the parameters alpha and r and time t:

<img src="https://latex.codecogs.com/svg.latex?p(t)&space;=&space;f(f_0,&space;r,&space;\alpha,&space;t)" title="p(t) = f(f_0, r, \alpha, t)" />

Parameters are shared across experiments (fixed effects), except for the initial conditions which might vary from experiment to experiment.

The measurement process is modelled with a normal distribution:
- In the case of an additive noise model:
<img src="https://latex.codecogs.com/svg.latex?p_\mathit{obs}(t)&space;\sim&space;\mathcal{N}&space;\big(&space;p(t),&space;\sigma^2&space;\big)" title="p_\mathit{obs}(t) \sim \mathcal{N} \big( p(t), \sigma^2 \big)" />

- In the case of an additive multiplicative noise model:
<img src="https://latex.codecogs.com/svg.latex?\log&space;\big(&space;p_\mathit{obs}(t)&space;\big)&space;\sim&space;\mathcal{N}&space;\Big(&space;\log&space;\big(&space;p(t)&space;\big),&space;\sigma^2&space;\Big)" title="\log \big( p_\mathit{obs}(t) \big) \sim \mathcal{N} \Big( \log \big( p(t) \big), \sigma^2 \Big)" />

### Priors

- The initial conditions of the different experiments are unknown and treated as parameters.
They are partially pooled with the hierarchical prior:

<img src="https://latex.codecogs.com/svg.latex?f_0&space;\sim&space;\mathcal{N}^&plus;(0,&space;\sigma_{f_0}^2)" title="f_0 \sim \mathcal{N}^+(0, \sigma_{f_0}^2)" />

- Competition rates are regularised with the hierarchical prior:

<img src="https://latex.codecogs.com/svg.latex?\alpha_{i,&space;j}&space;\sim&space;\mathcal{C}^&plus;(0,&space;\sigma_\alpha)" title="\alpha_{i, j} \sim \mathcal{C}^+(0, \sigma_\alpha)" />

Alternative choices could be considered, such as a Gaussian distribution (cf. ridge/L<sub>2</sub> regularisation), a Laplace distribution (cf. lasso/L<sub>1</sub> regularisation) or using a Horseshoe prior.

- Weakly informative priors are chosen for the remaining parameters.

## Results

Posterior predictive checks are conducted and the likelihood is stored to perform model selection with the loo package.
The additive and multiplicative models explain the data well but there is no strong evidence for choosing one model instead of the other (additive model is slightly favoured), although the experimental setup would suggest an additive noise.

Coefficients estimates for the two models (and those estimated in Gause's paper) are given below (error bars represents the 95% credible interval).
The original estimates are similar to those of the additive error model but slightly differ to the multiplicate error model.

![coefficient estimates](results.jpg)

## Future directions

- Regularisation alpha
- Sample k
- Random effects experiment level
- Alpha positive definite if systems stable
