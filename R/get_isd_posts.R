# Calculates the probability X>=x given values of lambda from the posterior (i.e., .epred). Allows for plotting the
# posterior of the isd

library(tidyverse)

get_isd_posts = function(data = NULL,
                         length = 300,
                         lambda = .epred){
  {{data}} %>%
    expand_grid(x.PLB = seq(unique(.$xmin), unique(.$xmax), length = 300)) %>%
    mutate(prob_yx = (1 - (x.PLB^({{lambda}} + 1) - (xmin^({{lambda}}+1)))/(xmax^({{lambda}} + 1) - (xmin^({{lambda}}+1)))))

}
