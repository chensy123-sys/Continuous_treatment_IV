# ------------------------------------------------------------
# Initialization
# ------------------------------------------------------------
rm(list = ls())                      # Clear workspace
source("cross_fitting.R")            # Load functions used for adaptive IV estimation

# ------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------
mydat <- read.csv("jpta_han.csv") %>% 
  select(-X) %>%                     # Remove index column
  mutate(A = pmax(pmin(A, 15), 8))   # Truncate education years to [8, 15]

table(mydat$A)                       # Check empirical distribution of A

# ------------------------------------------------------------
# Prepare storage for results
# ------------------------------------------------------------
A_value <- sort(unique(mydat$A))     # Support points of treatment A

res0 <- rep(0, length(A_value))      # SD of adaptive IV for L = 0 (Women)
res1 <- rep(0, length(A_value))      # SD of adaptive IV for L = 1 (Men)
res2 <- rep(0, length(A_value))      # SD of degenerate IV (baseline)

set.seed(2025)                       # Reproducibility

# ------------------------------------------------------------
# Compute variability of IV weights at each treatment value
# ------------------------------------------------------------
for (i in seq_along(A_value)) {
  
  # Adaptive RWF at target point A_value[i]
  delta <- delta_adap_func(mydat, target_point = A_value[i])
  
  # Adaptive RWF at target point A_value[i], where L is empty
  delta_deg <- delta_adap_func_deg(mydat, target_point = A_value[i])
  
  # Standard deviation of IV scores within sub-populations
  res0[i] <- sd(delta(mydat %>% filter(L == 0)))
  res1[i] <- sd(delta(mydat %>% filter(L == 1)))
  res2[i] <- sd(delta_deg(mydat))
  
  cat(i, " ")                        # Progress indicator
}

# ------------------------------------------------------------
# Prepare data for plotting
# ------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

df_plot <- data.frame(
  A_value = A_value,
  Total = res2,
  Women = res0,
  Men = res1
) %>%
  pivot_longer(
    cols = c(Total, Women, Men),
    names_to = "Group",
    values_to = "ChiSqDiv"
  )

# ------------------------------------------------------------
# Plot chi-square divergence across treatment levels
# ------------------------------------------------------------
p <- ggplot(df_plot, aes(x = A_value, y = ChiSqDiv, color = Group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = 0.05, 
    linetype = "dashed", 
    color = "gray40", 
    linewidth = 0.8
  ) +
  scale_color_manual(values = c(
    "Total" = "#1f77b4",
    "Women" = "#000000",
    "Men"   = "#d62728"
  )) +
  labs(
    x = "Years of education",
    y = expression(chi^2 ~ "divergence"),
    color = "Population"
  ) +
  coord_cartesian(ylim = c(0, 0.4)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y  = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text.x  = element_text(size = 18)
  )

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------
ggsave("data/Figure7.2(c).pdf", p, width = 6, height = 5)
print(p)
