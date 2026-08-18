# =============================================================================
# Project:     Submit DCC Jobs for WATT simulations with multiple treatments.
#              Uses methods from `simulation_methods_boot.R`
# =============================================================================
# Description: Description
#              Secondary Description Line
#
# Author:      Skyler Shapiro
# Updated:     2026-07-14
# =============================================================================
library(rslurm)
source("/hpc/home/sjs158/exp_bootvar/OWATT_functions.R")
# -----------------------------------------------------------------------------
# Tasks and settings
# -----------------------------------------------------------------------------
# my_tasks <- list(
#   trimming_05      = list(method = "trimming",         alpha = 0.05),
#   trimming_10      = list(method = "trimming",         alpha = 0.10),
#   trimming_15      = list(method = "trimming",         alpha = 0.15),
#   inf_trimming_04      = list(method = "influence_trimming", alpha = 0.04),
#   matching         = list(method = "matching_weights", alpha = 0) #,
# )
# 
# settings <- list(
#   good = list(kappas = c(0,0.08,0.2,0.1,0.13), N = 2000, alphas = c(.5,0.05,0.15,0.1, 0.1)),
#   adequate = list(kappas = c(0,.1,0.5,0.4,0.15), N = 2000, alphas = c(3,1.95,0.25,0.1, 1.8)),
#   poor = list(kappas = c(0,1,1.9,.87,0.8),N = 2000, alphas = c(5.5,1.05,0.15,0.2, 0.1))
# )

# -----------------------------------------------------------------------------
# Load results and compute metrics
# -----------------------------------------------------------------------------
source("/work/sjs158/multiple_WATT/0_params.R")
sjob   <- readRDS("sjob.rds")
res_all <- get_slurm_out(sjob, outtype = "raw", wrait = TRUE)
params  <- do.call(rbind, lapply(seq_along(settings), function(s)
  data.frame(seed = 1:1000, setting_idx = s)))

for (s in seq_along(settings)) {
  true_taus <- compute_true_tau(kappas = settings[[s]]$kappas,
                                alphas = settings[[s]]$alphas, tasks = my_tasks)
  res     <- res_all[which(params$setting_idx == s)]
  metrics_boot <- calculate_metrics(res, true_taus, se_source = "boot_se")
  saveRDS(metrics_boot, paste0("metrics_boot_s", s, ".rds"))
  cat("Setting:", names(settings)[s], "\n")
  print(metrics_boot)
}