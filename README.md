isdbayes: Bayesian hierarchical modeling of size spectra
================
Jeff Wesner

# Overview

This package allows the estimation of power law exponents using the
truncated (upper and lower) Pareto distribution (Wesner et al. 2023).
Specifically, it allows users to fit Bayesian (non)-linear hierarchical
models with a truncated Pareto likelihood using `brms` (Bürkner 2017).
The motivation for the package was to estimate power law exponents of
ecological size spectra using individual-level body size data in a
generalized mixed model framework. The likelihood for the truncated
Pareto used here was described in (Edwards et al. 2020). This package
translates that likelihood into `brms`.

# Installation

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

# Examples

``` r
# load these packages

library(dplyr)
library(tidyr)
library(here)
library(ggplot2)
library(tidybayes)
library(brms)
library(isdbayes)
```

## Fit individual samples

First, simulate some power law data using `rparetocounts()`. The code
below simulates 300 body sizes from a power law with exponent lambda =
-1.2, xmin = 1, and xmax = 1000.

``` r
# simulate data

dat = tibble(x = rparetocounts(n = 300,  lambda = -1.2,  xmin = 1, xmax = 1000)) |> 
  mutate(xmin = min(x),
         xmax = max(x),
         counts = 1)
```

The code above simulates data from a doubly-truncated Pareto and then
estimates xmin and xmax. It also adds a column for *counts.* If the data
all represent unique individual masses, then this column takes a value
of 1 for every body size. If the data have repeated sizes, then this
column can take an integer or double of the counts or densities of those
sizes. For example, data that are x = {1.9, 1.9, 1.8, 2.8, 2.8} could
either be analyzed with each body size assumed to be unique where counts
= {1, 1, 1, 1, 1} or it could be analyzed as x = {1.9, 1.8, 2.8} and
counts = {2, 1, 2}. The latter is a common format when there is a
density estimate associated with counts or a sampling effort.

Next estimate the power law exponent using `brms`. The model below
(`fit1`) is an intercept only model, where x are the body sizes and
counts, xmin, and xmax are included in `vreal()`. The use of `vreal` has
nothing to do with the model per se. It is simply required wording from
`brms` when including custom families. Similarly, `stanvars` is required
wording that contains the custom likelihood parameters. As long as
`isdbayes` is loaded, then `stanvars = stanvars` will work. It will stay
the same regardless of changes to the model structure (like new
predictors or varyaing intercepts).

``` r

fit1 = brm(x | vreal(counts, xmin, xmax) ~ 1, 
          data = dat,
          stanvars = stanvars,    # required for truncated Pareto
          family = paretocounts(),# required for truncated Pareto
          chains = 1, iter = 1000)
```

This example fits an intercept-only model to estimate the power-law
exponent. For more complex examples with fixed and hierarchical
predictors, see below.

## Simulate multiple size distributions

``` r

x1 = rparetocounts(lambda = -1.8) # `lambda` is required wording from brms. in this case it means the lambda exponent of the ISD
x2 = rparetocounts(lambda = -1.5)
x3 = rparetocounts(lambda = -1.2)

isd_data = tibble(x1 = x1,
                  x2 = x2,
                  x3 = x3) |> 
  pivot_longer(cols = everything(), names_to = "group", values_to = "x") |> 
  group_by(group) |> 
  mutate(xmin = min(x),
         xmax = max(x)) |> 
  group_by(group, x) |> 
  add_count(name = "counts")
```

## Fit multiple size distributions with a fixed factor

``` r
fit2 = brm(x | vreal(counts, xmin, xmax) ~ group, 
           data = isd_data,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000)
```

## Plot group posteriors

``` r
posts_group = fit2$data |> 
  distinct(group, xmin, xmax) |> 
  mutate(counts = 1) |> 
  add_epred_draws(fit2, re_formula = NA) 

posts_group |> 
  ggplot(aes(x = group, y = .epred)) + 
  stat_halfeye(scale = 0.2) + 
  geom_hline(yintercept = c(-1.8, -1.5, -1.2)) # known lambdas
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Fit multiple size distributions with a varying intercept

``` r
fit3 = brm(x | vreal(counts, xmin, xmax) ~ (1|group), 
           data = isd_data,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000)
```

## Plot varying intercepts

``` r
posts_varint = fit3$data |> 
  distinct(group, xmin, xmax) |> 
  mutate(counts = 1) |> 
  add_epred_draws(fit3, re_formula = NULL) 

posts_varint |> 
  ggplot(aes(x = group, y = .epred)) + 
  stat_halfeye(scale = 0.2) + 
  geom_hline(yintercept = c(-1.8, -1.5, -1.2)) # known lambdas
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Posterior predictive checks

After the model is fit, you can use built-in functions in brms to
perform model checking.

``` r
pp_check(fit2, type = "dens_overlay_grouped", group = "group") +
  scale_x_log10()
#> Using 10 posterior draws for ppc type 'dens_overlay_grouped' by default.
#> Warning in self$trans$transform(x): NaNs produced

#> Warning in self$trans$transform(x): NaNs produced
#> Warning: Removed 6 rows containing missing values (`geom_segment()`).
```

![](README_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-burkner2017brms" class="csl-entry">

Bürkner, P. C. 2017. “Brms: An r Package for Bayesian Multilevel Models
Using Stan.” *Journal of Statistical Software* 80: 1–28.

</div>

<div id="ref-edwards2020" class="csl-entry">

Edwards, A. M., J. P. W. Robinson, J. L. Blanchard, J. K. Baum, and M.
J. Plank. 2020. “Accounting for the Bin Structure of Data Removes Bias
When Fitting Size Spectra.” *Marine Ecology Progress Series* 636
(February): 19–33. <https://doi.org/10.3354/meps13230>.

</div>

<div id="ref-wesner2023bayesian" class="csl-entry">

Wesner, J. S, J. P. F. Pomeranz, J. R. Junker, and V. Gjoni. 2023.
“Bayesian Hierarchical Modeling of Size Spectra.” *bioRxiv*, 2023–02.

</div>

</div>
