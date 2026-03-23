library(np);library(dplyr)
env1 <- new.env()
env2 <- new.env()

load("data/LLKR_fit_5_10000.RData", envir = env1)
load("data/LLKR_fit_6_10000.RData", envir = env2)


LLKR_fit_list <- c(env1$LLKR_fit_list, env2$LLKR_fit_list)

Looptime <- length(LLKR_fit_list)
A_test<- seq(-0.25,0.25,0.01)
Y_AIPW <- matrix(0,Looptime, length(A_test))
Y_IPW <- matrix(0,Looptime, length(A_test))
Y_OR <- matrix(0,Looptime, length(A_test))

Y_AIPW_nuc <- matrix(0,Looptime, length(A_test))
Y_IPW_nuc <- matrix(0,Looptime, length(A_test))
Y_OR_nuc <- matrix(0,Looptime, length(A_test))



for (i in 1:Looptime) {
  Y_AIPW[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_AIPW,
    newdata = data.frame(treatment = A_test)
  )
  Y_IPW[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_IPW,
    newdata = data.frame(treatment = A_test)
  ) 
  Y_OR[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_OR,
    newdata = data.frame(treatment = A_test)
  )
  Y_AIPW_nuc[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_AIPW_nuc,
    newdata = data.frame(treatment = A_test)
  )
  Y_IPW_nuc[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_IPW_nuc,
    newdata = data.frame(treatment = A_test)
  ) 
  Y_OR_nuc[i,] <- predict(
    LLKR_fit_list[[i]]$LLKR_OR_nuc,
    newdata = data.frame(treatment = A_test)
  )
}



library(ggplot2);library(reshape2);library(grid);library(patchwork)



calc_mean_ci <- function(mat) {
  m <- apply(mat, 2, mean)
  se <- apply(mat, 2, sd)
  lower <- m - 1.96 * se
  upper <- m + 1.96 * se
  data.frame(mean = m, lower = lower, upper = upper)
}


df_list <- list(
  AIPW = calc_mean_ci(Y_AIPW),
  IPW = calc_mean_ci(Y_IPW),
  OR = calc_mean_ci(Y_OR),
  AIPW_nuc = calc_mean_ci(Y_AIPW_nuc),
  IPW_nuc = calc_mean_ci(Y_IPW_nuc),
  OR_nuc = calc_mean_ci(Y_OR_nuc),
  Truth = data.frame(mean = A_test, lower = A_test, upper = A_test)
)


Y_long <- do.call(rbind, lapply(names(df_list), function(method) {
  df <- df_list[[method]]
  data.frame(
    treatment = A_test,
    Method = method,
    mean = df$mean,
    lower = df$lower,
    upper = df$upper
  )
}))


Y_long_nuc <- subset(Y_long, Method %in% c("AIPW_nuc","IPW_nuc","OR_nuc","Truth"))

p_nuc <- ggplot(Y_long_nuc, aes(x = treatment, y = mean, color = Method, linetype = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c(
    "AIPW_nuc" = "blue",
    "IPW_nuc" = "red",
    "OR_nuc" = "green",
    "Truth" = "black"
  ),
  labels = c(
    "AIPW_nuc" = "AIPW NUC",
    "IPW_nuc" = "IPW NUC",
    "OR_nuc" = "OR NUC",
    "Truth" = "Truth",
    "AIPW" = "AIPW",
    "IPW" = "IPW",
    "OR" = "OR"
  )) +
  scale_fill_manual(values = c(
    "AIPW_nuc" = "blue",
    "IPW_nuc" = "red",
    "OR_nuc" = "green",
    "Truth" = "black"
  ),
  labels = c(
    "AIPW_nuc" = "AIPW NUC",
    "IPW_nuc" = "IPW NUC",
    "OR_nuc" = "OR NUC",
    "Truth" = "Truth",
    "AIPW" = "AIPW",
    "IPW" = "IPW",
    "OR" = "OR"
  )) +
  scale_linetype_manual(values = c(
    "AIPW_nuc" = "dashed",
    "IPW_nuc" = "dashed",
    "OR_nuc" = "dashed",
    "Truth" = "solid"
  ),
  labels = c(
    "AIPW_nuc" = "AIPW NUC",
    "IPW_nuc" = "IPW NUC",
    "OR_nuc" = "OR NUC",
    "Truth" = "Truth",
    "AIPW" = "AIPW",
    "IPW" = "IPW",
    "OR" = "OR"
  )) +
  labs(
    x = "Treatment",
    y = "Mean potential outcome curve",
    title = "NUC framework (n=10000)",
    color = "Method",
    linetype = "Method",
    fill = "Method"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) +
  ylim(-0.5, 0.5)





# 图2：No NUC + Truth
Y_long_no_nuc <- subset(Y_long, Method %in% c("AIPW","IPW","OR","Truth"))

p_no_nuc <- ggplot(Y_long_no_nuc, aes(x = treatment, y = mean, color = Method, linetype = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  scale_color_manual(values = c(
    "AIPW" = "blue",
    "IPW" = "red",
    "OR" = "green",
    "Truth" = "black"
  )) +
  scale_fill_manual(values = c(
    "AIPW" = "blue",
    "IPW" = "red",
    "OR" = "green",
    "Truth" = "black"
  )) +
  scale_linetype_manual(values = c(
    "AIPW" = "solid",
    "IPW" = "solid",
    "OR" = "solid",
    "Truth" = "solid"
  )) +
  labs(
    x = "Treatment",
    y = "Mean potential outcome curve",
    title = "IV framework (n=10000)",
    color = "Method",
    linetype = "Method",
    fill = "Method"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  ) +
  ylim(-0.5, 0.5)

ggsave("data/Figure6.2(f).pdf", plot = p_nuc + p_no_nuc + plot_layout(ncol = 2), width = 14, height = 5)
