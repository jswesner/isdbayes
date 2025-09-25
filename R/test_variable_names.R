stan_funs <- "
real paretocounts_lpdf(real Y, real lambda, real counts, real xmin, real xmax){
    if(lambda != -1)
    return(counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(Y)));
    else
    return(counts*-(log(log(xmax) - log(xmin)) + lambda*log(Y)));
}
"
stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

