
utils::globalVariables(c("x", "vreal2", "vreal3"))



#' Simulate data from a bounded power law
#' Simulates data from a bounded power law. This is basically a vectorized version of the rPLB() function
#' in Edwards et al. (2017) `sizeSpectra` package. The argument names are required to match the default
#' arguments for custom likelihoods in `brms`.
#'
#' @param n number of observations
#' @param mu vector of lambda (the power law exponent)
#' @param vreal2 xmin: the minimum body size of the sample or the minimum possible body size
#' @param vreal3 xmax: the maximum body size of the sample or the maximum possible body size
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' rparetocounts(n = 100, mu = -1.5, vreal2 = 1, vreal3 = 2000)
rparetocounts <- function(n = 300, mu = -1.2, vreal2 = 1, vreal3 = 1000) {
  samples <- numeric(n)
  {
    if(vreal2 <= 0 | vreal2 >= vreal3) stop("Parameters out of bounds in rPLB")
    u <- stats::runif(n)
    if(mu != -1){
     y <- ( u*vreal3^(mu+1) +  (1-u) * vreal2^(mu+1) ) ^ (1/(mu+1))
    } else
    { y <- vreal3^u * vreal2^(1-u)
    }
    return(y)
  }
  return(samples)
}


#' Bounded power law probability density function
#'
#' @param x body size value
#' @param mu vector of lambda (the power law exponent)
#' @param vreal2 xmin: the minimum body size of the sample or the minimum possible body size
#' @param vreal3 xmax: the maximum body size of the sample or the maximum possible body size
#'
#' @return a numeric vector of the value of the pdf given values of x, mu, xmin, and xmax.
#' @export
#'
#' @examples
#' xmin = 1
#' xmax = 1000
#' lambda = -1.5
#'
#' dparetocounts(x = 2, mu = lambda, vreal2 = xmin, vreal3 = xmax)
dparetocounts <- function(x, mu, vreal2, vreal3) {
  if (vreal2 <= 0 || vreal2 >= vreal3)
    stop("Parameters out of bounds in dPLB")

  if (x < vreal2 || x > vreal3)
    return(0)

  if (mu != -1) {
    density <- (mu + 1) * (x^(mu+1)) / (vreal3^(mu+1) - vreal2^(mu+1))
  } else {
    density <- x^(-2) / (vreal2 * log(vreal3/vreal2))
  }

  density
}


#' Generate custom log-likelihood for Stan
#'
#' This function is not called directly. It is only used to set up the log likelihood for
#' later use in `brms` via the family option within `brm()`
#'
#'
#' @param i see documentation from `get_dpar` within `brms`
#' @param prep see documentation from `get_dpar` within `brms`
#'
#' @return NA
#' @export
#'
#' @examples
#' NA
log_lik_paretocounts <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  vreal1 <- prep$data$vreal1[i]
  vreal2 <- prep$data$vreal2[i]
  vreal3 <- prep$data$vreal3[i]
  Y <- prep$data$Y[i]
  paretocounts_lpdf(Y, mu, vreal1, vreal2, vreal3)
}

#' Arrange data for posterior prediction
#'
#' This function is not called directly.
#'
#' @param i see documentation from `get_dpar` within `brms`
#' @param prep see documentation from `get_dpar` within `brms`
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
#' NA
posterior_predict_paretocounts <- function(i, prep, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  vreal2 = prep$data$vreal2[i]
  vreal3 = prep$data$vreal3[i]
  rparetocounts(prep$ndraws, mu, vreal2, vreal3)
}

#' Arrange data for posterior epred prediction
#'
#' This function is not called directly. It is only used to allow add_epred_draws() from the
#' `tidybayes` package to work with a `brmfit` object with a paretocounts_lpdf.
#'
#' @param prep see documentation from `get_dpar` within `brms`
#'
#' @return NA
#' @export
#'
#' @examples
#' NA
posterior_epred_paretocounts <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

stan_funs <- "
real paretocounts_lpdf(real Y, real mu, real vreal1, real vreal2, real vreal3){
    if(mu != -1)
    return(vreal1*(log((mu+1) / ( vreal3^(mu+1) - vreal2^(mu+1))) + mu*log(Y)));
    else
    return(vreal1*(log(log(vreal2) - log(vreal3)) + mu*log(Y)));
}
"
stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

#' Create custom paretocounts family
#'
#' This function is not called directly. It is only used to create the paretocounts() family option
#' in `brms`
#'
#'
#' @return NA
#' @export
#'
#' @examples
#' NA
paretocounts <- function(){custom_family(
  "paretocounts",
  dpars = c("mu"),
  links = c("identity"),
  lb = -Inf,
  ub = Inf,
  type = "real",
  vars = c("vreal1[n]",
           "vreal2[n]",
           "vreal3[n]")
)
}

#' Create paretocounts log posterior density function
#'
#' This function is not called directly. It is only used to create the paretocounts() family option
#' in `brms`
#'
#' @param Y vector of data (individual body sizes)
#' @param mu exponent of the bounded power law
#' @param vreal1 counts for each body size (numeric of integer)
#' @param vreal2 xmin: the minimum body size of the sample or the minimum possible body size
#' @param vreal3 xmax: the maximum body size of the sample or the maximum possible body size
#'
#' @return temporary function that is later vectorized using vectorize(paretocounts_lpdf_temp)
#' @export
#'
#' @examples
#' NA
paretocounts_lpdf_temp = function(Y, mu, vreal1, vreal2, vreal3){
  if(mu != -1){
    return(vreal1*(log((mu+1) / ( vreal3^(mu+1) - vreal2^(mu+1))) + mu*log(Y)));
  }
  if(mu == -1){
    return(vreal1*(log(log(vreal2) - log(vreal3)) + mu*log(Y)));
  }
}

paretocounts_lpdf = Vectorize(paretocounts_lpdf_temp)



