# =============================================================================
# Project:     Submit DCC Jobs for WATT simulations with multiple treatments.
#              Uses methods from `OWATT_functions.R`
# =============================================================================
#
# Author:      Skyler Shapiro
# Updated:     2026-07-22
# =============================================================================
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

source("/work/sjs158/multiple_WATT/0_params.R")

# Build one replicate row per setting [1000 individual jobs]
params <- do.call(rbind, lapply(seq_along(settings), function(s) {
  data.frame(seed = 1:1000, setting_idx = s)
}))

# Function each job will run
run_one <- function(seed, setting_idx) {
  source("/work/sjs158/multiple_WATT/OWATT_functions.R")
  s <- settings[[setting_idx]]
  print(s$kappas)
  single_sim(N = s$N, seed = seed, kappas = s$kappas, alphas = s$alphas,
             tasks = my_tasks, n_boot = 200)
}

sjob <- slurm_apply(
  f            = run_one,
  params       = params,
  jobname      = "watt_sim",
  nodes        = 2000,   
  cpus_per_node = 1,
  global_objects = c("settings", "my_tasks"),
  slurm_options = list(
    time       = "0:30:00",
    mem        = "2G",
    partition  = "biostat",
    account  = "biostat",
    `mail-user`  = "sjs158@duke.edu",
    `mail-type`  = "END,FAIL"
  ),
  submit = TRUE
)

saveRDS(sjob, "sjob.rds")