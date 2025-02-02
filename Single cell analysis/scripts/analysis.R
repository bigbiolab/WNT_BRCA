# Single Cell Analysis
# Author: Muntasim Fuad

# Load required packages
library(tidyverse)
library(openxlsx)

# Load datasets 
wnt2 <- read.csv("Single cell analysis/data/TISCH_WNT2_heatmap.csv")
wnt7b <- read.csv("Single cell analysis/data/TISCH_WNT7B_heatmap.csv")
wnt11 <- read.csv("Single cell analysis/data/TISCH_WNT11_heatmap.csv")

# Filter expression based values less than zero
wnt2 <- wnt2 |> filter(Values > 0)
wnt7b <- wnt7b |> filter(Values > 0)
wnt11 <- wnt11 |> filter(Values > 0)

# ----------------------------------WNT2----------------------------------------
# Generate bubble heatmaps
wnt2_plot <- ggplot(wnt2, aes(x = Datasets, 
                              y = Category, 
                              size = Values, 
                              fill = Values)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient(low = "#fc9272", high = "#e34a33") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    panel.border = element_rect(color = "#636363", fill = NA, linewidth = 1),
    legend.position = "right"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  labs(title = "WNT2",
       x = "",
       y = "",
       size = "log(TPM/10+1)",
       fill = "log(TPM/10+1)")

# Save the plots
ggsave(filename = "Single cell analysis/figures/WNT2.png",
       plot = wnt2_plot,
       width = 4,
       height = 8,
       units = "in",
       dpi = 300,
       bg = "white")

# ----------------------------------WNT7B---------------------------------------
# Generate bubble heatmaps
wnt7b_plot <- ggplot(wnt7b, aes(x = Datasets, 
                              y = Category, 
                              size = Values, 
                              fill = Values)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient(low = "#fc9272", high = "#e34a33") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    panel.border = element_rect(color = "#636363", fill = NA, linewidth = 1),
    legend.position = "right")+
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
      labs(title = "WNT7B",
       x = "",
       y = "",
       size = "log(TPM/10+1)",
       fill = "log(TPM/10+1)")

# Save the plots
ggsave(filename = "Single cell analysis/figures/WNT7B.png",
       plot = wnt7b_plot,
       width = 4,
       height = 8,
       units = "in",
       dpi = 300,
       bg = "white")

# ----------------------------------WNT11---------------------------------------
# Generate bubble heatmaps
wnt11_plot <- ggplot(wnt11, aes(x = Datasets, 
                                y = Category, 
                                size = Values, 
                                fill = Values)) +
  geom_point(shape = 21, color = "#636363") +
  scale_size(range = c(2, 10)) +
  scale_fill_gradient(low = "#fc9272", high = "#de2d26") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               vjust = 1),
    panel.border = element_rect(color = "#636363", fill = NA, linewidth = 1),
    legend.position = "right"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)) +
  labs(title = "WNT11",
       x = "",
       y = "",
       size = "log(TPM/10+1)",
       fill = "log(TPM/10+1)")

# Save the plots
ggsave(filename = "Single cell analysis/figures/WNT11.png",
       plot = wnt11_plot,
       width = 4,
       height = 8,
       units = "in",
       dpi = 300,
       bg = "white")
# ------------------------------------------------------------------------------
# Export outputs
write.xlsx(wnt2,"Single cell analysis/outputs/WNT2.xlsx")
write.xlsx(wnt7b,"Single cell analysis/outputs/WNT7B.xlsx")
write.xlsx(wnt11,"Single cell analysis/outputs/WNT11.xlsx")