library(rslurm)
# Specify models
# my_tasks <- list(
#   trimming_05      = list(method = "trimming",         alpha = 0.05),
#   trimming_10      = list(method = "trimming",         alpha = 0.10),
#   trimming_15      = list(method = "trimming",         alpha = 0.15),
#   inf_trimming_04      = list(method = "influence_trimming", alpha = 0.04),
#   matching         = list(method = "matching_weights", alpha = 0)
# )
# 
# # Specify overlap settings
# settings <- list(
#   poor = list(kappas =   c(0,0.4,0.1,0.4,.1), N = 2000, alphas = c(3,.95,.25, 2.5, .3)),
#   good =  list(kappas =  c(0,.2,.3, .2,.3), N = 2000, alphas =  c(2.2,.95,.25, 1, .3))
# )
source("/work/sjs158/multiple_watt/0_params.R")
truth_params <- expand.grid(rep = 1:80, setting_idx = seq_along(settings))


run_truth_rep <- function(rep, setting_idx) {
  source("/hpc/home/sjs158/exp_bootvar/OWATT_functions.R")
 N <- 100000 
s <- settings[[setting_idx]]
  set.seed(42 + rep)
  X  <- make_data(N)
  ps <- make_treatment(X, N, kappas = s$kappas, alphas = s$alphas, report_ps = TRUE)
  Y  <- make_outcome(X, ps[, 1], N)
  lapply(my_tasks, function(t)
    sapply(2:5, function(i)
      WATT(Y[, 6], ps[, 1], X, reference_group = 1, treatment_group = i,
           method = t$method, alpha = t$alpha, true_e = ps[, -1])))
}

truth_job <- slurm_apply(run_truth_rep, truth_params, jobname = "watt_truth",
                         nodes = nrow(truth_params), cpus_per_node = 1,
                         global_objects = c("settings", "my_tasks"),
                         slurm_options = list(
                           time       = "0:30:00",
                           mem        = "2G",
                           partition  = "biostat",
                           account  = "biostat",
                           `mail-user`  = "sjs158@duke.edu",
                           `mail-type`  = "END,FAIL"
                         ), submit = TRUE)

saveRDS(list(truth_job = truth_job, truth_params = truth_params,
             my_tasks = my_tasks, settings = settings),
        "job_meta.rds")
