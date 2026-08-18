my_tasks <- list(
  trimming_05      = list(method = "trimming",         alpha = 0.05),
  trimming_10      = list(method = "trimming",         alpha = 0.10),
  trimming_15      = list(method = "trimming",         alpha = 0.15),
  inf_trimming_04      = list(method = "influence_trimming", alpha = 0.04),
  matching         = list(method = "matching_weights", alpha = 0)
)

# Specify overlap settings
settings <- list(
  poor = list(kappas =   c(0,0.4,0.1,0.4,.1), N = 2000, alphas = c(3,.95,.25, 2.5, .3)),
  good =  list(kappas =  c(0,.2,.3, .2,.3), N = 2000, alphas =  c(2.2,.95,.25, 1, .3))
)
