library(np);library(dplyr);rm(list=ls()); source('cross_fitting2.R')
library(tidyr);library(ggplot2)
env9 <- new.env(); load("data/JTPA_9.RDATA", envir = env9)
env11 <- new.env(); load("data/JTPA_11.RDATA", envir = env11)
env12 <- new.env(); load("data/JTPA_12.RDATA", envir = env12)
env15 <- new.env(); load("data/JTPA_15.RDATA", envir = env15)
env9_deg <- new.env(); load("data/JTPA_9_deg.RDATA", envir = env9_deg)
env11_deg <- new.env(); load("data/JTPA_11_deg.RDATA", envir = env11_deg)
env12_deg <- new.env(); load("data/JTPA_12_deg.RDATA", envir = env12_deg)
env15_deg <- new.env(); load("data/JTPA_15_deg.RDATA", envir = env15_deg)

# Number of Monte Carlo repetitions
n_mc <- length(env9$fit_CF_list)

# Treatment grid for smooth curves
A_grid <- seq(8, 15, 0.1)

# Treatment values used for spline fitting
# treat_subset <- c(8, 9, 10, 13, 14, 15)
treat_subset <- c(8, 9, 10,12, 13, 14, 15); alpha <- 0.00

############################################################
# Helper function: smoothing spline prediction
############################################################

smooth_predict <- function(data, y_var, A_grid, treat_subset) {
  
  subdata <- data %>%
    filter(treatment %in% treat_subset) 
  
  fit <- smooth.spline(
    x = subdata$treatment,
    y = subdata[[y_var]],tol = 1e-6
  )
  
  predict(fit, x = A_grid)$y
}



############################################################
# Storage containers
############################################################

score_summary_list <- vector("list", n_mc)
curve_nuc <- matrix(NA, n_mc, length(A_grid))
curve_nuc_deg <- matrix(NA, n_mc, length(A_grid))
curve_IV <- matrix(NA, n_mc, length(A_grid))
curve_IV_deg <- matrix(NA, n_mc, length(A_grid))


############################################################
# Monte Carlo loop
############################################################

for (i in seq_len(n_mc)) {
  
  score_9 <- env9$fit_CF_list[[i]]$score
  score_11 <- env11$fit_CF_list[[i]]$score
  score_12 <- env12$fit_CF_list[[i]]$score
  score_15 <- env15$fit_CF_list[[i]]$score
  score_merge <- rbind(
    score_9 %>% filter(treatment %in% c(8,9,10)),
    score_11 %>% filter(treatment %in% c(11)),
    score_12 %>% filter(treatment %in% c(12)),
    score_15 %>% filter(treatment %in% c(13,14,15))
  )
  
  score_9_deg <- env9_deg$fit_CF_list[[i]]$score
  score_11_deg <- env11_deg$fit_CF_list[[i]]$score
  score_12_deg <- env12_deg$fit_CF_list[[i]]$score
  score_15_deg <- env15_deg$fit_CF_list[[i]]$score
  score_merge_deg <- rbind(
    score_9_deg %>% filter(treatment %in% c(8,9,10)),
    score_11_deg %>% filter(treatment %in% c(11)),
    score_12_deg %>% filter(treatment %in% c(12)),
    score_15_deg %>% filter(treatment %in% c(13,14,15))
  )
  
  # Smoothed curves
  curve_nuc[i, ] <- smooth_predict(
    score_9, "score_AIPW_nuc", A_grid, treat_subset
  )
  curve_nuc_deg[i, ] <- smooth_predict(
    score_9_deg, "score_AIPW_nuc", A_grid, treat_subset
  )
  curve_IV[i, ] <- smooth_predict(
    score_merge %>%
      filter(
        score_AIPW > quantile(score_AIPW, alpha),
        score_AIPW < quantile(score_AIPW, 1-alpha)
      ), "score_AIPW", A_grid, treat_subset
  )
  curve_IV_deg[i, ] <- smooth_predict(
    score_merge_deg %>%
      filter(
        score_AIPW > quantile(score_AIPW, alpha),
        score_AIPW < quantile(score_AIPW, 1-alpha)
      ), "score_AIPW", A_grid, treat_subset
  )
  
  if (i %% 100 == 0) cat(i, " ")
}

############################################################
# Construct plotting dataset
############################################################

data_plot <- data.frame(
  treatment = A_grid,
  
  est_NUC = colMeans(curve_nuc),
  lower_NUC = apply(curve_nuc, 2, quantile, 0.025),
  upper_NUC = apply(curve_nuc, 2, quantile, 0.975),
  
  est_IV = colMeans(curve_IV),
  lower_IV = apply(curve_IV, 2, quantile, 0.025),
  upper_IV = apply(curve_IV, 2, quantile, 0.975),
  
  est_IV_deg = colMeans(curve_IV_deg),
  lower_IV_deg = apply(curve_IV_deg, 2, quantile, 0.025),
  upper_IV_deg = apply(curve_IV_deg, 2, quantile, 0.975)
)

############################################################
# Convert to long format for ggplot
############################################################


data_long <- data_plot %>%
  pivot_longer(
    cols = -treatment,
    names_to = c(".value", "group"),
    names_pattern = "(est|lower|upper)_(NUC|IV_deg|IV)"
  )

############################################################
# Plot estimated ADRF curves
############################################################

p <- ggplot(data_long,
            aes(x = treatment, y = est,
                color = group, fill = group)) +
  
  geom_line(linewidth = 1) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    color = NA
  ) +
  
  labs(
    # title = "",
    x = "Years of education",
    y = "Estimated prev-program earning",
    color = "Estimator:",
    fill = "Estimator:"
  ) +
  
  scale_color_manual(
    values = c(
      "NUC"     = "#7570b3",
      "IV"      = "#d95f02",
      "IV_deg"  = "#e7298a"
    ),
    labels = c(
      "NUC"     = "NUC (L: sex)",
      "IV"      = "IV (L: sex)",
      "IV_deg"  = "IV (L: empty)"
    )
  ) +
  scale_fill_manual(
    values = c(
      "NUC"     = "#7570b3",
      "IV"      = "#d95f02",
      "IV_deg"  = "#e7298a"
    ),
    labels = c(
      "NUC"     = "NUC (L: sex)",
      "IV"      = "IV (L: sex)",
      "IV_deg"  = "IV (L: empty)"
    )
  )+
  
  theme_minimal(base_size = 12) +
  
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 10),     
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 16)  
  ) + guides(
    color = guide_legend(nrow = 2),
    fill  = guide_legend(nrow = 2)
  )
print(p)
ggsave('data/Figure7.1(c).pdf',p,width = 6,height = 5)





# library(np);library(dplyr);rm(list=ls()); source('cross_fitting2.R')
# library(tidyr);library(ggplot2)
# env9 <- new.env(); load("data/JTPA_9.RDATA", envir = env9)
# env11 <- new.env(); load("data/JTPA_11.RDATA", envir = env11)
# env12 <- new.env(); load("data/JTPA_12.RDATA", envir = env12)
# env15 <- new.env(); load("data/JTPA_15.RDATA", envir = env15)
# env9_deg <- new.env(); load("data/JTPA_9_deg.RDATA", envir = env9_deg)
# env11_deg <- new.env(); load("data/JTPA_11_deg.RDATA", envir = env11_deg)
# env12_deg <- new.env(); load("data/JTPA_12_deg.RDATA", envir = env12_deg)
# env15_deg <- new.env(); load("data/JTPA_15_deg.RDATA", envir = env15_deg)
# 
# # Number of Monte Carlo repetitions
# n_mc <- length(env9$fit_CF_list)
# 
# # Treatment grid for smooth curves
# A_grid <- seq(8, 15, 0.1)
# 
# # Treatment values used for spline fitting
# # treat_subset <- c(8, 9, 10, 13, 14, 15)
# treat_subset <- c(8, 9, 10,12, 13, 14, 15); alpha <- 0.00
# 
# ############################################################
# # Helper function: smoothing spline prediction
# ############################################################
# 
# smooth_predict <- function(data, y_var, A_grid, treat_subset) {
#   
#   subdata <- data %>%
#     filter(treatment %in% treat_subset) 
#   
#   fit <- smooth.spline(
#     x = subdata$treatment,
#     y = subdata[[y_var]],tol = 1e-6
#   )
#   
#   predict(fit, x = A_grid)$y
# }
# 
# 
# 
# ############################################################
# # Storage containers
# ############################################################
# 
# score_summary_list <- vector("list", n_mc)
# curve_nuc <- matrix(NA, n_mc, length(A_grid))
# curve_nuc_deg <- matrix(NA, n_mc, length(A_grid))
# curve_IV <- matrix(NA, n_mc, length(A_grid))
# curve_IV_deg <- matrix(NA, n_mc, length(A_grid))
# 
# 
# ############################################################
# # Monte Carlo loop
# ############################################################
# 
# for (i in seq_len(n_mc)) {
#   
#   score_9 <- env9$fit_CF_list[[i]]$score
#   score_11 <- env11$fit_CF_list[[i]]$score
#   score_12 <- env12$fit_CF_list[[i]]$score
#   score_15 <- env15$fit_CF_list[[i]]$score
#   score_merge <- rbind(
#     score_9 %>% filter(treatment %in% c(8,9,10)),
#     score_11 %>% filter(treatment %in% c(11)),
#     score_12 %>% filter(treatment %in% c(12)),
#     score_15 %>% filter(treatment %in% c(13,14,15))
#   )
#   
#   score_9_deg <- env9_deg$fit_CF_list[[i]]$score
#   score_11_deg <- env11_deg$fit_CF_list[[i]]$score
#   score_12_deg <- env12_deg$fit_CF_list[[i]]$score
#   score_15_deg <- env15_deg$fit_CF_list[[i]]$score
#   score_merge_deg <- rbind(
#     score_9_deg %>% filter(treatment %in% c(8,9,10)),
#     score_11_deg %>% filter(treatment %in% c(11)),
#     score_12_deg %>% filter(treatment %in% c(12)),
#     score_15_deg %>% filter(treatment %in% c(13,14,15))
#   )
#   
#   # Smoothed curves
#   curve_nuc[i, ] <- smooth_predict(
#     score_9, "score_AIPW_nuc", A_grid, treat_subset
#   )
#   curve_nuc_deg[i, ] <- smooth_predict(
#     score_9_deg, "score_AIPW_nuc", A_grid, treat_subset
#   )
#   curve_IV[i, ] <- smooth_predict(
#     score_merge %>%
#       filter(
#         score_AIPW > quantile(score_AIPW, alpha),
#         score_AIPW < quantile(score_AIPW, 1-alpha)
#       ), "score_AIPW", A_grid, treat_subset
#   )
#   curve_IV_deg[i, ] <- smooth_predict(
#     score_merge_deg %>%
#       filter(
#         score_AIPW > quantile(score_AIPW, alpha),
#         score_AIPW < quantile(score_AIPW, 1-alpha)
#       ), "score_AIPW", A_grid, treat_subset
#   )
#   
#   if (i %% 100 == 0) cat(i, " ")
# }
# 
# ############################################################
# # Construct plotting dataset
# ############################################################
# 
# data_plot <- data.frame(
#   treatment = A_grid,
#   
#   est_NUC = colMeans(curve_nuc),
#   lower_NUC = apply(curve_nuc, 2, quantile, 0.025),
#   upper_NUC = apply(curve_nuc, 2, quantile, 0.975),
#   
#   est_NUC_deg = colMeans(curve_nuc),
#   lower_NUC_deg = apply(curve_nuc_deg, 2, quantile, 0.025),
#   upper_NUC_deg = apply(curve_nuc_deg, 2, quantile, 0.975),
#   
#   est_IV = colMeans(curve_IV),
#   lower_IV = apply(curve_IV, 2, quantile, 0.025),
#   upper_IV = apply(curve_IV, 2, quantile, 0.975),
#   
#   est_IV_deg = colMeans(curve_IV_deg),
#   lower_IV_deg = apply(curve_IV_deg, 2, quantile, 0.025),
#   upper_IV_deg = apply(curve_IV_deg, 2, quantile, 0.975)
# )
# 
# ############################################################
# # Convert to long format for ggplot
# ############################################################
# 
# 
# data_long <- data_plot %>%
#   pivot_longer(
#     cols = -treatment,
#     names_to = c(".value", "group"),
#     names_pattern = "(est|lower|upper)_(NUC_deg|NUC|IV_deg|IV)"
#   )
# 
# ############################################################
# # Plot estimated ADRF curves
# ############################################################
# 
# p <- ggplot(data_long,
#             aes(x = treatment, y = est,
#                 color = group, fill = group)) +
#   
#   geom_line(linewidth = 1) +
#   
#   geom_ribbon(
#     aes(ymin = lower, ymax = upper),
#     alpha = 0.2,
#     color = NA
#   ) +
#   
#   labs(
#     # title = "",
#     x = "Years of education",
#     y = "Estimated prev-program earning",
#     color = "Estimator:",
#     fill = "Estimator:"
#   ) +
#   
#   scale_color_manual(
#     values = c(
#       "NUC_deg" = "#1b9e77",
#       "NUC"     = "#7570b3",
#       "IV"      = "#d95f02",
#       "IV_deg"  = "#e7298a"
#     ),
#     labels = c(
#       "NUC_deg" = "NUC (L: empty)",
#       "NUC"     = "NUC (L: sex)",
#       "IV"      = "IV (L: sex)",
#       "IV_deg"  = "IV (L: empty)"
#     )
#   ) +
#   scale_fill_manual(
#     values = c(
#       "NUC_deg" = "#1b9e77",
#       "NUC"     = "#7570b3",
#       "IV"      = "#d95f02",
#       "IV_deg"  = "#e7298a"
#     ),
#     labels = c(
#       "NUC_deg" = "NUC (L: empty)",
#       "NUC"     = "NUC (L: sex)",
#       "IV"      = "IV (L: sex)",
#       "IV_deg"  = "IV (L: empty)"
#     )
#   )+
#   
#   theme_minimal(base_size = 12) +
#   
#   theme(
#     legend.position = "top",
#     panel.grid.minor = element_blank(),
#     axis.title = element_text(size = 18),
#     axis.text = element_text(size = 18),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(size = 10),     
#     axis.title.x = element_text(size = 16),   
#     axis.title.y = element_text(size = 16),
#     legend.title = element_text(size = 16),   
#     legend.text = element_text(size = 16)  
#   ) + guides(
#     color = guide_legend(nrow = 2),
#     fill  = guide_legend(nrow = 2)
#   )
# print(p)
# ggsave('data/Figure7.1(c).pdf',p,width = 6,height = 5)
# 
# 
