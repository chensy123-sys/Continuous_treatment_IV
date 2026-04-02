library(ggplot2)    # Load ggplot2 for plotting
library(dplyr)      # Load dplyr for data manipulation
library(latex2exp)  # Load latex2exp to use LaTeX expressions in plots

# Define the sequence of x values from 0 to 1 with step 0.005
x <- seq(0, 1, 0.005)

# Prepare the data frame for plotting
plot_data <- data.frame(
  x = rep(x, 2),  # Repeat x values for each group
  y = c(
    dbeta(x, 2, 6),  # Beta distribution for Z=0
    dbeta(x, 6, 2)   # Beta distribution for Z=1
  ),
  group = factor(
    rep(c("Z=0", "Z=1"), each = length(x)),  # Create group labels
    levels = c("Z=0", "Z=1")                 # Set factor levels explicitly
  )
)

# Create ggplot object
p <- ggplot(plot_data, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1.4) +  # Draw lines with specified width
  scale_color_manual(
    values = c("steelblue", "firebrick"),                        # Custom colors
    labels = TeX(c("$p_{A|Z}(a|0)$", "$p_{A|Z}(a|1)$")),        # LaTeX labels for legend
    name = TeX("$p_{A|Z}(a|Z):$")                                # Legend title in LaTeX
  ) +
  labs(
    title = "",             # No main title
    x = "A (Treatment)",    # X-axis label
    y = TeX("Density")      # Y-axis label using LaTeX
  ) +
  theme_minimal(base_size = 15) +  # Minimal theme with base font size
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),  # Center title
    axis.title.x = element_text(size = 16),            # X-axis title size
    axis.title.y = element_text(size = 16),            # Y-axis title size
    axis.text = element_text(size = 16),              # Axis text size
    legend.position = "top",                           # Position legend at top
    legend.title = element_text(size = 14),           # Legend title font size
    legend.text = element_text(size = 13)             # Legend text font size
  )

# Save the plot as a PDF with width=4 inches and height=4 inches
ggsave('data/Figure2.1(a).pdf', p, width = 4, height = 4)