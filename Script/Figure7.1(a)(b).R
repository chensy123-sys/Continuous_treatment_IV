library(dplyr);library(np);library(mgcv);library(coin);library(ggplot2)
rm(list = ls()); source('cross_fitting2.R'); h1 <- 1; h2 <- 2

mydat <- read.csv("jpta_han.csv") %>% select(-X) %>%
  mutate(
    # Y = as.numeric(Y>median(Y)),
    A = pmax(pmin(A,15),8)
  )

A_test <- seq(8,15)
# A_test <- c(9,14)


# L = sex
plot_data <- data.frame()
for (a in A_test) {
  delta <- delta_adap_func(mydat, target_point = a,h = h1, constriant = NULL)
  res <- URWF_test(data = mydat, func = delta, region = seq(8,15), constriant = 'L', h=h2,divide = 2)
  temp <- data.frame(A_test = a, region = res$region, pval = res$pvals)
  plot_data <- bind_rows(plot_data, temp)
  print(a)
}

p1 <- ggplot(plot_data, aes(x = region, y = log(pval), color = factor(paste0('a=',A_test)))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "black") +
  labs(x = "Years of education",y = "logarithm of p-value",color = expression(paste(hat(pi)[a], " for "))) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 10),     
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 16)   
  )

# L = empty
plot_data <- data.frame()
for (a in A_test) {
  delta <- delta_adap_func(mydat, target_point = a,h = h1, constriant = NULL)
  res <- URWF_test(data = mydat, func = delta, region = seq(8,15), h=h2,constriant = NULL)
  temp <- data.frame(A_test = a, region = res$region, pval = res$pvals)
  plot_data <- bind_rows(plot_data, temp)
  print(a)
}

p2 <- ggplot(plot_data, aes(x = region, y = log(pval), color = factor(paste0('a=',A_test)))) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = log(0.05), linetype = "dashed", color = "black") +
  labs(x = "Years of education",y = "logarithm of p-value",color = expression(paste(hat(pi)[a], " for "))) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 10),     
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 16),   
    legend.text = element_text(size = 16)   
  )
print(p1+p2)
ggsave('data/Figure7.1(a).pdf',p1,device = 'pdf',width = 6,height = 5)
ggsave('data/Figure7.1(b).pdf',p2,device = 'pdf',width = 6,height = 5)

