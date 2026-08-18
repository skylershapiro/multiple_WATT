library(writexl)
source("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/sandvar/OWATT_functions.R")

# sandwich_1 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/sandvar/metrics_boot_s1.rds")
# sandwich_2 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/sandvar/metrics_boot_s2.rds")
# sandwich_3 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/sandvar/metrics_boot_s3.rds")

sandwich_results <- list("good" = sandwich_1, "adequate" = sandwich_2, "poor" = sandwich_3)
# 
boot_1 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/bootvar/0701/metrics_boot_s1.rds")
boot_2 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/bootvar/0701/metrics_boot_s2.rds")
boot_3 <- readRDS("/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/bootvar/0701/metrics_boot_s3.rds")

boot_results <- list("good" = boot_1, "adequate" = boot_2, "poor" = boot_3)

write_xlsx(sandwich_results, path = "/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/sandvar/sandwich_variance_results.xlsx")
write_xlsx(boot_results, path = "/Users/skylershapiro/duke/DCRI/OWATT_simulations_skyler/simulation_results/bootvar/bootstrap_variance_results.xlsx")


# >>>>>>>>>>>>>> GENERATE OVERLAP PLOTS <<<<<<<<<<<<<<
# good_kappas <- c(0, 0.2, 0.2, 0.2, 0.2)
# adequate_kappas <- c(0, 0.8, 0.8, 0.8, 0.8)
# poor_kappas <-  c(0, 1.5, 1.5, 1.5, 1.5)
#plot_overlap(good_kappas)

library(ggplot2)
library(RColorBrewer)

plot_method_comparison <- function(df, metric = c("ARBias_pct", "RMSE", "Avg_Estimated_SE")) {
  
  metric <- match.arg(metric)
  
  ylab <- switch(metric,
                 ARBias_pct = "Absolute Relative Bias (%)",
                 RMSE       = "Root Mean Squared Error (RMSE)",
                 Avg_Estimated_SE   = "Average Estimated Standard Error")
  
  trim_colors  <- brewer.pal(4, "Blues")[2:4]
  trunc_colors <- brewer.pal(4, "Purples")[2:4]
  
  method_colors <- c(
    "conventional_att" = "#2F4F4F",
    "overlap"          = "#E69F00",
    "inf_trimming_04"  = "#A87C4F",
    "matching"         = "#CC3399",
    
    "trimming_05"   = trim_colors[1],
    "trimming_10"   = trim_colors[2],
    "trimming_15"   = "#4E7CA8",
    
    "truncation_05" = trunc_colors[1],
    "truncation_10" = trunc_colors[2],
    "truncation_15" = "#74689E"
  )
  
  ggplot(df, aes(x = Treatment, y = .data[[metric]], fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8),
             width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = method_colors) +
    scale_x_discrete(labels = c("ATT_2_vs_1" = expression(lambda[12]),
                                "ATT_3_vs_1" = expression(lambda[13]),
                                "ATT_4_vs_1" = expression(lambda[14]),
                                "ATT_5_vs_1" = expression(lambda[15]))) +
    theme_minimal(base_family = "serif") +
    labs(x = "Pairwise Contrast", y = ylab) +
    theme(
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dashed")
    )
}

plot_method_comparison_facet <- function(df_list, metric = "ARBias_pct",
                                         facet_labels = NULL) {
  if (is.null(facet_labels)) facet_labels <- names(df_list)
  if (is.null(facet_labels)) facet_labels <- paste0("Sandwich ", seq_along(df_list))
  ylab <- switch(metric,
                 ARBias_pct = "Absolute Relative Bias (%)",
                 RMSE       = "Root Mean Squared Error (RMSE)",
                 Avg_Estimated_SE   = "Average Estimated Standard Error",
                 CP = "Coverage Probability")
  
  combined <- bind_rows(
    lapply(seq_along(df_list), function(i) {
      d <- df_list[[i]]
      d$Scenario <- facet_labels[i]
      d
    })
  )
  combined$Scenario <- factor(combined$Scenario, levels = facet_labels)
  trim_colors  <- brewer.pal(4, "Blues")[2:4]
  trunc_colors <- brewer.pal(4, "Purples")[2:4]
  method_colors <- c(
    "conventional_att" = "#2F4F4F",
    "overlap"          = "#E69F00",
    "inf_trimming_04"  = "#A87C4F",
    "matching"         = "#CC3399",
    "trimming_05"   = trim_colors[1],
    "trimming_10"   = trim_colors[2],
    "trimming_15"   = "#4E7CA8",
    "truncation_05" = trunc_colors[1],
    "truncation_10" = trunc_colors[2],
    "truncation_15" = "#74689E"
  )
  if (metric == "CP"){ggplot(combined, aes(x = Treatment, y = .data[[metric]], color = Method)) +
      geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = 1.00), 
                fill = "grey90", color = NA) +
      geom_point(size = 3, position = position_dodge(width = 0.6)) +
      scale_color_manual(values = method_colors) +
      geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", linewidth = 0.5) +
      scale_x_discrete(labels = c("ATT_2_vs_1" = expression(lambda[12]),
                                  "ATT_3_vs_1" = expression(lambda[13]),
                                  "ATT_4_vs_1" = expression(lambda[14]),
                                  "ATT_5_vs_1" = expression(lambda[15]))) +
      facet_wrap(~ Scenario, nrow = 1) +
      theme_minimal(base_family = "serif") +
      labs(x = "Pairwise Contrast", y = ylab) +
      theme(
        axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(color = "black", size = 12, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        strip.text = element_text(color = "black", size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dashed"))}
  else{
  ggplot(combined, aes(x = Treatment, y = .data[[metric]], fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(0.8),
             width = 0.9, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = method_colors) +
    scale_x_discrete(labels = c("ATT_2_vs_1" = expression(lambda[12]),
                                "ATT_3_vs_1" = expression(lambda[13]),
                                "ATT_4_vs_1" = expression(lambda[14]),
                                "ATT_5_vs_1" = expression(lambda[15]))) +
    facet_wrap(~ Scenario, nrow = 1) +
    theme_minimal(base_family = "serif") +
    labs(x = "Pairwise Contrast", y = ylab) +
    theme(
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 12, face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      strip.text = element_text(color = "black", size = 12, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4, linetype = "dashed")
    )}
}

# Parameters
alpha <- 0.2
k <- 1000 
e_j <- seq(0, 1, length.out = 500)

# Sigmoid function definition
sigmoid <- function(x, k) { 1 / (1 + exp(-k * x)) }

# Approximation: 
# The difference of two sigmoids creates the "box" shape
h_approx <- sigmoid(e_j - alpha, k) - sigmoid(e_j - (1 - alpha), k)

# Plotting
plot(e_j, h_approx, type = "l", lwd = 3, col = "red",
     main = "Sigmoid approximation of trimming selection function (k=1000)",
     xlab = expression(e[j]), ylab = "h (approx)")

# Overlaying original function for comparison
lines(e_j, (e_j >= alpha) & (e_j <= (1 - alpha)), col = "blue", lty = 2, lwd = 2)

legend("topright", legend=c("Sigmoid Approx", "Original"), 
       col=c("red", "blue"), lty=c(1, 2), lwd=c(3, 2))

# Usage:
plot_method_comparison(boot_1, metric = "ARBias_pct")
# plot_method_comparison(sandwich_1, metric = "RMSE")



plot_method_comparison_facet(
  list(boot_1, boot_2, boot_3),
  metric = "ARBias_pct",
  facet_labels = c("Good Overlap", "Adequate Overlap", "Poor Overlap")
)


plot_method_comparison_facet(
  list(boot_1, boot_2, boot_3),
  metric = "RMSE",
  facet_labels = c("Good Overlap", "Adequate Overlap", "Poor Overlap")
)


plot_method_comparison_facet(
  list(boot_1, boot_2, boot_3),
  metric = "Avg_Estimated_SE",
  facet_labels = c("Good Overlap", "Adequate Overlap", "Poor Overlap")
)

plot_method_comparison_facet(
  list(boot_1, boot_2, boot_3),
  metric = "CP",
  facet_labels = c("Good Overlap", "Adequate Overlap", "Poor Overlap")
)
