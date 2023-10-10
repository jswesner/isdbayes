isdbayes: Bayesian hierarchical modeling of size spectra
================
Jeff Wesner
2023-10-04

## Overview

This package allows the estimation of power law exponents using the
truncated (upper and lower) Pareto distribution. Specifically, it allows
users to fit Bayesian (non)-linear hierarchical models with a truncated
Pareto likelihood using `brms` (Bürkner 2017). The motivation for the
package was to estimate power law exponents of ecological size spectra
using individual-level body size data in a generalized mixed model
framework. The likelihood for the truncated Pareto used here was
described in (Edwards et al. 2020). This package translates that
likelihood into `brms`.

## Installation

This package requires installation of `brms` and `rstan`, which itself
requires installation of a C++ toolchain.

1)  Go to <https://mc-stan.org/users/interfaces/rstan.html> and follow
    the instructions to install `rstan` and configure the C++ toolchain.

2)  Install the latest version of `brms` with install.packages(“brms”).

3)  Install `isdbayes` from github using `devtools`:

``` r
# requires an installation of devtools

devtools::install_github("jswesner/isdbayes")
```

## Example

`isdbayes` amends the stanvars and family options in brms models to
accept the truncated Pareto. First, simulate some power law data using
the `rparetocounts()` function. The code below simulates 300 body sizes
from a power law with exponent (mu) = -1.2, xmin (vreal2) = 1, and xmax
(vreal3) = 1000. The options are called “mu”, “vreal2”, and “vreal3”
instead of “lambda”, “xmin”, and “xmax” to fit with the generic naming
requirements of `brms` for custom family distributions.

``` r
library(isdbayes)

# simulate data
dat = tibble(x = rparetocounts(n = 300,  mu = -1.2,  vreal2 = 1, vreal3 = 1000)) %>% 
  mutate(xmin = min(x),
         xmax = max(x),
         counts = 1)
```

The code above simulates data from a doubly-truncated Pareto and then
estimates xmin and xmax. It also adds a column for *counts.* If the data
all represent unique individual masses, then this column takes a value
of 1 for every body size. If the data have repeated sizes, then this
column can take an integer of the counts of those sizes. For example,
data that are x = {1.9, 1.9, 1.8, 2.8, 2.8} could either be analyzed
with each body size assumed to be unique where counts = {1, 1, 1, 1, 1}
or it could be analyzed as x = {1.9, 1.8, 2.8} and counts = {2, 1, 2}.
The latter is a common format when there is a density estimate
associated with counts or a sampling effort.

Next estimate the power law exponent using `brms`.

``` r
library(brms)

m1 = brm(x | vreal(counts, xmin, xmax) ~ 1, 
          data = dat,
          stanvars = stanvars,    # required for truncated Pareto
          family = paretocounts(),# required for truncated Pareto
          chains = 1, iter = 1000)
```

This example fits an intercept-only model to estimate the power-law
exponent. For more complex examples with fixed and hierarchical
predictors, [see the
vignette](https://github.com/jswesner/isdbayes/blob/master/vignettes/isdbayes_vignette.Rmd).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-burkner2017brms" class="csl-entry">

Bürkner, Paul-Christian. 2017. “Brms: An r Package for Bayesian
Multilevel Models Using Stan.” *Journal of Statistical Software* 80:
1–28.

</div>

<div id="ref-edwards2020" class="csl-entry">

Edwards, Am, Jpw Robinson, Jl Blanchard, Jk Baum, and Mj Plank. 2020.
“Accounting for the Bin Structure of Data Removes Bias When Fitting Size
Spectra.” *Marine Ecology Progress Series* 636 (February): 19–33.
<https://doi.org/10.3354/meps13230>.

</div>

</div>
