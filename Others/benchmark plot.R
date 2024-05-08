# Load necessary libraries
library(ggplot2)
library(scales)
library(extrafont)

# Load the reshaped data from a CSV
plot_data <- read.csv("Benchmark_Allocation_Rates_reshaped.csv")

# Register the fonts with R
loadfonts(device = "win")

# Plot using ggplot2 with "Computer Modern" font and enlarged axis values
plot <- ggplot(plot_data, aes(x = Age, y = Benchmark_Rate, color = Fund_Type, linetype = Fund_Type)) +
  geom_line() +
  scale_x_continuous(breaks = seq(25, 65, by = 5), limits = c(25, 65)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.1),
    labels = percent_format(scale = 100)
  ) +
  labs(x = "Age", y = "Allocation to Risky Asset", color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "CMU Serif"),  # Set the font family to "Computer Modern"
    axis.text = element_text(size = 12.5, family = "CMU Serif"),  # Enlarge the axis values and set the font
    legend.text = element_text(size = 12.5, family = "CMU Serif")  # Enlarge the axis values and set the font
  )

print(plot)



