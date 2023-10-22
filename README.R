knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## # requires an installation of devtools
##
## devtools::install_github("jswesner/isdbayes")
##

## library(isdbayes)
##
## # simulate data
 dat = tibble(x = rparetocounts(n = 300,  mu = -1.2,  vreal2 = 1, vreal3 = 1000)) %>%
   mutate(xmin = min(x),
          xmax = max(x),
          counts = 1)


library(brms)

fit1 = brm(x | vreal(counts, xmin, xmax) ~ 1,
           data = dat,
           stanvars = stanvars,    # required for truncated Pareto
           family = paretocounts(),# required for truncated Pareto
           chains = 1, iter = 1000)
saveRDS(fit1, file = "models/fit1.rds")

library(isdbayes)
library(tidyverse)
library(brms)

x1 = rparetocounts(mu = -1.8) # `mu` is required wording from brms. in this case it means the lambda exponent of the ISD
x2 = rparetocounts(mu = -1.5)
x3 = rparetocounts(mu = -1.2)

isd_data = tibble(x1 = x1,
                  x2 = x2,
                  x3 = x3) %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "x") %>%
  group_by(group) %>%
  mutate(xmin = min(x),
         xmax = max(x)) %>%
  group_by(group, x) %>%
  add_count(name = "counts")

fit2 = brm(x | vreal(counts, xmin, xmax) ~ group,
           data = isd_data,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000,
           file_refit = "on_change",
           file = "models/fit2.rds")

posts_group = fit2$data %>%
  distinct(group, xmin, xmax) %>%
  mutate(counts = 1) %>%
  tidybayes::add_epred_draws(fit2, re_formula = NA)

posts_group %>%
  ggplot(aes(x = group, y = .epred)) +
  tidybayes::stat_halfeye(scale = 0.2) +
  geom_hline(yintercept = c(-1.8, -1.5, -1.2)) # known lambdas

fit3 = update(fit2, formula = . ~ (1|group))

saveRDS(fit3, file = "models/fit3.rds")

posts_varint = fit3$data %>%
  distinct(group, xmin, xmax) %>%
  mutate(counts = 1) %>%
  tidybayes::add_epred_draws(fit3, re_formula = NULL)

posts_varint %>%
  ggplot(aes(x = group, y = .epred)) +
  tidybayes::stat_halfeye(scale = 0.2) +
  geom_hline(yintercept = c(-1.8, -1.5, -1.2)) # known lambdas

pp_check(fit2, type = "dens_overlay_grouped", group = "group")

WAIC(fit2)
WAIC(fit3)

