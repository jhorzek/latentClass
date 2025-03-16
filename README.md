
<!-- README.md is generated from README.Rmd. Please edit that file -->

# latentClass

<!-- badges: start -->

<!-- badges: end -->

latentClass is an R package that provides a vanilla implementation of
latent class models with continuous variables (Gaussian) and categorical
variables.

latentClass tries to identify subgroups of individuals within the data
set. The central assumption is that, once we know the class an
individual is part of, all items are independent (local independence).
The likelihood of the model is given by:

$$L = \prod_{i=1}^N \prod_{c=1}^C \pi_cp(x_{i1}|c)p(x_{i2}|c)\cdots p(x_{iP}|c),$$

where $N$ is the sample size, $C$ is the number of classes, $\pi_c$ is
the probability of class c, $p(x_{i1}|c)$ is the likelihood of observing
$x_{i1}$ given that person $i$ is in class $c$, and $P$ is the number of
items.

Currently, latentClass supports normal (Gaussian) items and categorical
items:

- Gaussian items: The likelihood for normal (Gaussian) items is given by
  $p(x) = \frac{1}{\sqrt{2\pi\sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$
  (see <https://en.wikipedia.org/wiki/Normal_distribution> for more
  details).
- Categorical items: The likelihood for categorical items is given by
  $p(x) = p_1^{[x=1]}p_2^{[x=2]}\cdots p_k^{[x=k]}$ (see
  <https://en.m.wikipedia.org/wiki/Categorical_distribution> for more
  details).

Parameters are estimated with an Expectation Maximization optimizer.
Additionally, latentClass supports sample weights (see Murphy & Scrucca,
2012).

latentClass is one of multiple packages that implement latent class
models in R. The implementation of categorical items is heavily inspired
by the poLCA package (Drew et al., 2011). Similarly, the implementation
of Gaussian items is inspired by mclust (Scrucca, et al., 2023). The
implementation of the Expectation Maximization optimizer follows that
found in Blume (2002). Overall latentClass is a very vanilla
implementation of latent class models, lacking many of the more advanced
features of other packages (e.g., no latent class regressions).

## Installation

You can install the development version of latentClass from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jhorzek/latentClass")
```

## Example

``` r
library(latentClass)
data <- simulate_latent_class_data()
model <- latentClass(data = data,
                    # Define the number of classes:
                    n_classes = 3,
                    # Specify which of the items are categorical:
                    categorical = categorical(items = c("cat_1", "cat_2")),
                    # Specify which of the items are Gaussian:
                    normal = normal(items = c("norm_1", "norm_2")))

summary(model)
#> #### Latent Class Model Results #####
#> 
#> Model settings:
#> -------------- 
#> - Number of classes:     3
#> - Categorical variables: cat_1, cat_2
#> - Normal variables:      norm_1, norm_2
#> 
#> Estimation:
#> ---------- 
#> - Model converged: Yes 
#> - Estimation time: 0.003 seconds 
#> 
#> Fit measures:
#> ------------ 
#> - Prameters:       60
#> - Observations:    500
#> - log-Likelihood: -2552.877
#> - BIC:             5478.63
#> - AIC:             5225.753
#> 
#> Estimates:
#> --------- 
#> cat_1:
#>     class_1   class_2    class_3
#> 1 0.1888938 0.3003217 0.80443911
#> 2 0.6329966 0.2953494 0.10533002
#> 3 0.1781097 0.4043288 0.09023086
#> 
#> cat_2:
#>     class_1   class_2   class_3
#> 1 0.1435871 0.6085802 0.3009979
#> 2 0.8564129 0.3914198 0.6990021
#> 
#> norm_1:
#>         class_1     class_2    class_3
#> mean -0.1485329 -0.05203985 0.06275182
#> sd    0.9819375  0.98193750 0.98193750
#> 
#> norm_2:
#>         class_1    class_2   class_3
#> mean 0.03731104 -2.9170354 2.9679840
#> sd   0.85786628  0.8578663 0.8578663
```

## References:

- Blume, M. (2002). Expectation maximization: A gentle introduction.
  Technical University of Munich Institute for Computer Science.
- Drew A. Linzer, Jeffrey B. Lewis (2011). poLCA: An R Package for
  Polytomous Variable Latent Class Analysis. Journal of Statistical
  Software, 42(10), 1-29. URL <https://www.jstatsoft.org/v42/i10/>
- Murphy, T. B., & Scrucca, L. (2012). Using Weights in mclust.
- Scrucca L, Fraley C, Murphy TB, Raftery AE (2023). *Model-Based
  Clustering, Classification, and Density Estimation Using mclust in R*.
  Chapman and Hall/CRC. ISBN 978-1032234953, <doi:10.1201/9781003277965>
  <https://doi.org/10.1201/9781003277965>,
  <https://mclust-org.github.io/book/>
