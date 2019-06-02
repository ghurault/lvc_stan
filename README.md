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

![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bdp%7D%7Bdt%7D%20%3D%20p%20%28r%20-%20%5Calpha%20p%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

![equation](http://www.sciweavers.org/tex2img.php?eq=p%28t%29%20%3D%20%5Cfrac%7Bk%7D%7B1%20%2B%20C%20e%5E%7B-%20r%20t%7D%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

- r is the growth parameter
- alpha is the (self) competition rate
- k = r / alpha is the carrying capacity (steady state)
- C is a constant determined by initial conditions

For two species:

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Bcases%7D%0A%5Cfrac%7Bdp_1%7D%7Bdt%7D%20%3D%20p_1%20%28r_1%20-%20%5Calpha_%7B1%2C%201%7D%20p_1%20-%20%5Calpha_%7B1%2C2%7D%20p_2%29%20%5C%5C%0A%5Cfrac%7Bdp_2%7D%7Bdt%7D%20%3D%20p_2%20%28r_2%20-%20%5Calpha_%7B2%2C%201%7D%20p_1%20-%20%5Calpha_%7B2%2C%202%7D%20p_2%29%0A%5Cend%7Bcases%7D&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\begin{cases}\frac{dp_1}{dt} = p_1 (r_1 - \alpha_{1, 1} p_1 - \alpha_{1,2} p_2) \\\frac{dp_2}{dt} = p_2 (r_2 - \alpha_{2, 1} p_1 - \alpha_{2, 2} p_2)\end{cases}" width="260" height="53" />

More generally, for N species (`*` is the element-wise multiplication):

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%7Bd%5Cboldsymbol%7Bp%7D%7D%7Bdt%7D%20%3D%20%5Cboldsymbol%7Bp%7D%20%2A%20%28%5Cboldsymbol%7Br%7D%20-%20A%20%5Cboldsymbol%7Bp%7D%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\frac{d\boldsymbol{p}}{dt} = \boldsymbol{p} * (\boldsymbol{r} - A \boldsymbol{p})" width="150" height="43" />

### Statistical model

The solutions of the ODEs can expressed as function of the initial conditions f_0, the parameters alpha and r and time t:

<img src="http://www.sciweavers.org/tex2img.php?eq=p%28t%29%20%3D%20f%28p_0%2C%20r%2C%20%5Calpha%2C%20t%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="p(t) = f(f_0, r, \alpha, t)" width="150" height="19" />

The measurement processed is modelled with a normal distribution:
- In the case of an additive noise model:
<img src="http://www.sciweavers.org/tex2img.php?eq=p_%5Cmathit%7Bobs%7D%28t%29%20%5Csim%20%5Cmathcal%7BN%7D%20%5Cbig%28%20p%28t%29%2C%20%5Csigma%5E2%20%5Cbig%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="p_\mathit{obs}(t) \sim \mathcal{N} \big( p(t), \sigma^2 \big)" width="175" height="22" />
- In the case of an additive multiplicative noise model:
<img src="http://www.sciweavers.org/tex2img.php?eq=%5Clog%20%5Cbig%28%20p_%5Cmathit%7Bobs%7D%28t%29%20%5Cbig%29%20%5Csim%20%5Cmathcal%7BN%7D%20%5CBig%28%20%5Clog%20%5Cbig%28%20p%28t%29%20%5Cbig%29%2C%20%5Csigma%5E2%20%5CBig%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\log \big( p_\mathit{obs}(t) \big) \sim \mathcal{N} \Big( \log \big( p(t) \big), \sigma^2 \Big)" width="272" height="32" />

### Priors

The initial conditions of the different experiments are partially pooled with the hierarchical prior:

![equation](http://www.sciweavers.org/tex2img.php?eq=f_0%20%5Csim%20%5Cmathcal%7BN%7D%5E%2B%280%2C%20%5Csigma_%7Bf_0%7D%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

Competition rates are regularised with the hierarchical prior:

<img src="http://www.sciweavers.org/tex2img.php?eq=%5Calpha_%7Bi%2C%20j%7D%20%5Csim%20%5Cmathcal%7BC%7D%5E%2B%280%2C%20%5Csigma_%5Calpha%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="\alpha_{i, j} \sim \mathcal{C}^+(0, \sigma_\alpha)" width="122" height="22" />

Alternative choices could be considered, such as a Gaussian distribution (cf. ridge/L<sub>2</sub> regularisation), a Laplace distribution (cf. lasso/L<sub>1</sub> regularisation) or using a Horseshoe prior.

Weakly informative priors were chosen for the remaining parameters.

## Results

Posterior predictive checks are conducted and the likelihood is stored to perform model selection with the loo package.

![coefficient estimates](results.jpg)
