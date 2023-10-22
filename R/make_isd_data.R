make_isd_data = function(model = NULL, group = NULL, x = x){
  {{model}}$data %>%
    group_by(across({{group}})) %>%
    arrange(-{{x}}) %>%
    mutate(n_yx = 1:max(row_number()))
}
