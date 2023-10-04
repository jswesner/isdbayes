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
described in (Edwards et al. 2020). We translated it to `rstan` (Stan
Development Team, n.d.) and `brms`, allowing users to use the helpful
tools in those packages.

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

<div id="ref-rstan" class="csl-entry">

Stan Development Team. n.d. “RStan: The R Interface to Stan.”
<https://mc-stan.org/>.

</div>

</div>
