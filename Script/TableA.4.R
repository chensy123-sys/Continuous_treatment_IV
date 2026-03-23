library(np);library(dplyr);rm(list=ls()); source('cross_fitting2.R')
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
# Storage containers
############################################################

score_summary_list <- vector("list", n_mc)
score_summary_list_deg <- vector("list", n_mc)



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
  
  # Summary statistics by treatment
  score_summary_list[[i]] <- cbind(
    score_9 %>%
      group_by(treatment) %>%
      summarise(
        AIPW_9 = mean(score_AIPW),
        sd_AIPW_9 = sd(score_AIPW) / sqrt(n()),
        AIPW_nuc = mean(score_AIPW_nuc),
        sd_AIPW_nuc = sd(score_AIPW_nuc) / sqrt(n()),
        .groups = "drop"
      ),
    score_11 %>%
      group_by(treatment) %>%
      summarise(
        AIPW_11 = mean(score_AIPW),
        sd_AIPW_11 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_12 %>%
      group_by(treatment) %>%
      summarise(
        AIPW_12 = mean(score_AIPW),
        sd_AIPW_12 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_15 %>%
      group_by(treatment) %>%
      summarise(
        AIPW_15 = mean(score_AIPW),
        sd_AIPW_15 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_merge %>%
      group_by(treatment) %>%
      summarise(
        AIPW_merge = mean(score_AIPW),
        sd_AIPW_merge = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment) 
  )
  
  
  # Summary statistics by treatment
  score_summary_list_deg[[i]] <- cbind(
    score_9_deg %>%
      group_by(treatment) %>%
      summarise(
        AIPW_9 = mean(score_AIPW),
        sd_AIPW_9 = sd(score_AIPW) / sqrt(n()),
        AIPW_nuc = mean(score_AIPW_nuc),
        sd_AIPW_nuc = sd(score_AIPW_nuc) / sqrt(n()),
        .groups = "drop"
      ),
    score_11_deg %>%
      group_by(treatment) %>%
      summarise(
        AIPW_11 = mean(score_AIPW),
        sd_AIPW_11 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_12_deg %>%
      group_by(treatment) %>%
      summarise(
        AIPW_12 = mean(score_AIPW),
        sd_AIPW_12 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_15_deg %>%
      group_by(treatment) %>%
      summarise(
        AIPW_15 = mean(score_AIPW),
        sd_AIPW_15 = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment),
    score_merge_deg %>%
      group_by(treatment) %>%
      summarise(
        AIPW_merge = mean(score_AIPW),
        sd_AIPW_merge = sd(score_AIPW) / sqrt(n()),
        .groups = "drop"
      ) %>% select(-treatment) 
  )
  
  if (i %% 100 == 0) cat(i, " ")
}


############################################################
# Combine Monte Carlo summaries
############################################################


score_summary <- bind_rows(score_summary_list)
score_summary_deg <- bind_rows(score_summary_list_deg)

summary_table <- score_summary %>%
  group_by(treatment) %>%
  summarise(
    est_9 = mean(AIPW_9),
    sd_9 = sd(AIPW_9),
    est_12 = mean(AIPW_12),
    sd_12 = sd(AIPW_12),
    est_15 = mean(AIPW_15),
    sd_15 = sd(AIPW_15),
    est_NUC = mean(AIPW_nuc),
    sd_NUC = sd(AIPW_nuc),
    est_IV = mean(AIPW_merge),
    sd_IV = sd(AIPW_merge),
    .groups = "drop"
  ) %>% filter(treatment %in% treat_subset)

summary_table_deg <- score_summary_deg %>%
  group_by(treatment) %>%
  summarise(
    est_9 = mean(AIPW_9),
    sd_9 = sd(AIPW_9),
    est_12 = mean(AIPW_12),
    sd_12 = sd(AIPW_12),
    est_15 = mean(AIPW_15),
    sd_15 = sd(AIPW_15),
    est_NUC = mean(AIPW_nuc),
    sd_NUC = sd(AIPW_nuc),
    est_IV = mean(AIPW_merge),
    sd_IV = sd(AIPW_merge),
    .groups = "drop"
  ) %>% filter(treatment %in% treat_subset)

res <- data.frame(t(rbind(summary_table_deg,summary_table)))[-1,]
colnames(res) <- c(treat_subset,treat_subset)
write.csv(round(res), file = "data/TableA.4.csv")
# library(xtable)
# xtable(round(res),)