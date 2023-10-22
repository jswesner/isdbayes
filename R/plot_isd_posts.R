# Calculates the probability X>=x given values of lambda from the posterior (i.e., .epred).
# Then plots the ISD with log10 axes.

plot_isd_posts = function(model = NULL, group = NULL, ndraws = 500,
                         nvalues = 100, re_formula = NA, countname = "counts",
                         xmin = "xmin", xmax = "xmax"){

  posts = model$data %>%
    distinct(across({{group}}), xmin, xmax) %>%
    mutate({{countname}} := 1) %>%
    add_epred_draws(model, re_formula = re_formula, ndraws = ndraws)

  plot_data = model$data %>%
    group_by(across({{group}}), xmin, xmax) %>%
    add_tally() %>%
    distinct(across({{group}}), n) %>%
    group_by(across({{group}}), n) %>%
    mutate(min_x_group = min(xmin), max_x_group = max(xmax)) %>%
    group_modify(function(df, .keep = "unused") {
      # Create a new column "x" with a logarithmic sequence between min_x_group and max_x_group
      df <- df %>%
        expand(x = 10^seq(log10(min_x_group), log10(max_x_group), length.out = nvalues))

      return(df)
    }) %>%
    left_join(posts,
              relationship = "many-to-many") %>%
    mutate(prob_yx = (1 - (x^(.epred + 1) - (xmin^(.epred+1)))/((xmax)^(.epred + 1) - (xmin^(.epred+1)))),
           n_yx = prob_yx*n)

  plot_data %>%
    group_by(across({{group}}), x) %>%
    median_qi(n_yx) %>%
    ggplot(aes(x = x, y = n_yx, group = {{group}})) +
    geom_line() +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.4) +
    scale_y_log10() +
    scale_x_log10() +
    coord_cartesian(ylim = c(1, NA))

}
