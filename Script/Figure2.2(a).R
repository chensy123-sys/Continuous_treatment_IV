library(dplyr)      # Data manipulation
library(np)         # Nonparametric estimation
library(mgcv)       # Generalized Additive Models (GAM)
library(fields)     # Tools for spatial data and visualization
library(latex2exp)  # Use LaTeX expressions in plots
library(ggplot2)    # Plotting

rm(list = ls())     # Clear the workspace

# Function to simulate data
sim_data <- function(n){
  mydat <- data.frame(L = rnorm(n, sd = 1)) %>%  # Covariate L
    mutate(
      U = -0.5 * L + rnorm(n, sd = 2),           # Unobserved confounder
      Z = -0.5 * L + rnorm(n, sd = 2),           # Instrument variable
      sel = rbinom(n, 1, 0.5),                   # Selector for treatment
      A = ifelse(sel, -0.5 * L + Z, -0.5 * L + U) + rnorm(n), # Treatment
    ) %>%
    select(-sel)                         # Remove unused columns
  return(mydat)
}

n <- 5e4
set.seed(2025)
mydat <- sim_data(n)  # Generate simulated data

# Fit GAM models
rho_model <- gam(Z ~ s(L), data = mydat, method = "REML")          # E[Z|L] model
rho <- as.vector(predict(rho_model, newdata = mydat, type = 'response'))

kappa_model <- gam(Z ~ te(A, L), data = mydat, method = "REML")   # E[Z|A,L] model
kappa <- as.vector(predict(kappa_model, newdata = mydat %>% mutate(A = -0.3), type='response')) - rho

# Compute 5% lower quantile of |kappa|
lower_kappa <- quantile(abs(kappa), 0.05)

# Construct grid for plotting
new_L <- seq(-1, 1, 0.01)        # Grid for L
i_seq <- seq(-3, 3, 0.5)         # Grid for A

# Predict E[Z|L] on the L grid
rho <- as.vector(predict(rho_model, newdata = data.frame(L = new_L), type = 'response'))

# Expand grid for plotting
plot_data <- expand.grid(L = new_L, A = i_seq)
plot_data$rho <- rep(rho, times = length(i_seq))                     # Add rho
plot_data$kappa <- as.vector(predict(kappa_model, newdata = plot_data, type = 'response')) - plot_data$rho  # Compute kappa

# Plot kappa vs L for different A
p <- ggplot(plot_data, aes(x = L, y = kappa, group = A, color = A)) +
  geom_line(linewidth = 0.75) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7) +  # Add horizontal reference line at y=0
  scale_color_gradientn(colors = colorRampPalette(c("blue", "red"))(200)) +  # Color gradient for A
  labs(
    title = "",
    x = "L",
    y = TeX("$E[Z| A=a,L] - E[Z| L]$"),
    color = "A"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 16)
  )

# Save the plot as a PDF
ggsave('data/Figure2.2(a).pdf', p, width = 6, height = 5)