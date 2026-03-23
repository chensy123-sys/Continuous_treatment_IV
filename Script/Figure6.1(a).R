library(dplyr);library(np);library(mgcv);library(coin);library(ggplot2);rm(list = ls());
source('cross_fitting2.R')
set.seed(2025)
delta <- delta_adap_func(sim_data(1e4,theta_func = theta), target_point=0.5)
n <- 10000; mydat <- sim_data(n,theta_func = theta)
# A_test <- seq(-0.95, 0.95, length.out = 5)
A_test <- c(-0.5,0,0.5); region = seq(-0.75, 0.75, 0.01)
plot_data <- data.frame()
for (a in A_test) {
  delta <- delta_adap_func(mydat, target_point = a)
  res <- URWF_test(data = mydat, func = delta, region = region, constriant = 'L',divide = 1)
  temp <- data.frame(A_test = a, region = res$region, pval = res$pvals)
  plot_data <- bind_rows(plot_data, temp)
}

p <- ggplot(plot_data, aes(x = region, y = log(pval), color = factor(paste0('a=',A_test)))) +
  geom_line() +
  geom_point(size = 1) +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Treatment region: [-0.75, 0.75]",
    y = "Logarithm of p-value",
    color = expression(paste(hat(pi)[a], " for "))
  ) +
  theme_minimal() + scale_x_continuous(breaks = seq(-0.75, 0.75, by = 0.25)) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 10),     
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 16)     
  )
print(p)
ggsave(filename = 'data/Figure6.1(a).pdf',plot = p,device = 'pdf',width = 7,height = 5)






