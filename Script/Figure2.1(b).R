library(ggplot2)  # Load ggplot2 for plotting
library(dplyr)    # Load dplyr for data manipulation

# Define the sequence of x values from 0 to 1 with step 0.005
x <- seq(0, 1, 0.005)

# Prepare the data frame for plotting
plot_data <- data.frame(
  x = rep(x, 3),  # Repeat x for each group
  y = c(
    dbeta(x, 2, 6),  # Beta distribution for Z=0
    dbeta(x, 4, 4),  # Beta distribution for Z=1
    dbeta(x, 6, 2)   # Beta distribution for Z=2
  ),
  group = factor(
    rep(c("Z=0", "Z=1", "Z=2"), each = length(x)),  # Create group labels
    levels = c("Z=0", "Z=1", "Z=2")                 # Set factor levels
  )
)

# Create ggplot object
p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1.4) +  # Plot lines with specified width
  scale_color_manual(
    values = c("steelblue", "darkorange", "firebrick"),  # Set custom colors
    labels = TeX(c("$p_{A|Z}(a|0)$", "$p_{A|Z}(a|1)$", "$p_{A|Z}(a|2)$")),  # LaTeX labels for legend
    name = TeX("$p_{A|Z}(a|Z):$")  # Legend title in LaTeX
  ) +
  labs(
    title = "",          # No main title
    x = "A (Treatment)", # X-axis label
    y = TeX("Density")   # Y-axis label with LaTeX
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  # Centered title
    axis.title.x = element_text(size = 16),             # X-axis title size
    axis.title.y = element_text(size = 16),             # Y-axis title size
    axis.text = element_text(size = 16),               # Axis text size
    legend.position = "top",                            # Position legend at top
    legend.title = element_text(size = 14),            # Legend title size
    legend.text = element_text(size = 13)              # Legend text size
  )

# Save the plot as PDF with width=4 inches and height=4 inches
ggsave('data/Figure2.1(b).pdf', p, width = 4, height = 4)
