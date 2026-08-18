library(rslurm)
source("/work/sjs158/multiple_WATT/OWATT_functions.R")  # for calculate_metrics only
source("/work/sjs158/multiple_WATT/0_params.R")
meta <- readRDS("job_meta.rds")
sjob <- readRDS("sjob.rds")
list2env(meta, envir = environment())

sim_job      <- sjob$sim_job

# Build one replicate row per setting [1000 individual jobs]
sim_params <- do.call(rbind, lapply(seq_along(settings), function(s) {
  data.frame(seed = 1:1000, setting_idx = s)
}))

truth_job    <- meta$truth_job
truth_params <- meta$truth_params
my_tasks     <- meta$my_tasks
settings     <- meta$settings

sim_res   <- get_slurm_out(sjob, outtype = "raw", wait = TRUE)
truth_res <- get_slurm_out(truth_job, outtype = "raw", wait = TRUE)

for (s in seq_along(settings)) {
  
  true_taus <- lapply(names(my_tasks), function(nm) {
    reps <- truth_res[truth_params$setting_idx == s]
    vals <- do.call(rbind, lapply(reps, `[[`, nm))
    setNames(colMeans(vals), paste0("ATT_", 2:5, "_vs_1"))
  })
  names(true_taus) <- names(my_tasks)
  
  res <- sim_res[sim_params$setting_idx == s]
  metrics_boot <- calculate_metrics(res, true_taus, se_source = "boot_se")
  
  saveRDS(metrics_boot, paste0("metrics_boot_s", s, ".rds"))
  cat("Setting:", names(settings)[s], "\n")
  print(metrics_boot)
}
