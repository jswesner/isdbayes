library(dplyr)

# make data for isd plot with x = body sizes and y = prob X>=x or num X>=x
get_isd_data = function(data = NULL, resp = Y){
  data %>% arrange(desc({{resp}})) %>%
    mutate(rownumber = row_number(),
           num_yx = 1:max(rownumber),
           prob_yx = seq(0, 1, length.out = max(rownumber)))
}
