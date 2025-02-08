# Breast Cancer Cell Lines & Gene Effect Score Analysis
# Author: Muntasim Fuad

# Load required pacakges
library(tidyverse)
library(cowplot)

# Load gene effect data
wnt2 <- read.csv("Cancer Cell line data analysis//data/WNT2.csv")
wnt7b <- read.csv("Cancer Cell line data analysis//data/WNT7B.csv")
wnt11 <- read.csv("Cancer Cell line data analysis//data/WNT11.csv")


# .........................WNT2.......................................
# Categorize the gene effect scores as Negative or Positive
wnt2$Gene.effect <- ifelse(wnt2$Gene.effect.score < 0, 
                                     "negative", "positive")

# Create bar plot
wnt2_plot <- ggplot(wnt2, aes(x = Category, 
                                y = Gene.effect.score,
                                fill = Gene.effect)) +
  geom_bar(stat = "identity", width=0.2) + 
  theme_minimal() +
  labs(title = "Gene Effect Score for WNT2 in Breast Cancer Cell Lines",
       x = "Breast Cancer Cell Lines",
       y = "Gene Effect Score") +
  scale_fill_manual(values = c("negative" = "#9A133DFF", "positive" = "#1A318BFF")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 14, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   colour = "#636363" ),
        axis.text.y = element_text(size = 14,
                                   colour = "#00010DFF" ), 
        plot.title = element_text(size = 20, face = "bold", color = "#071140FF"), 
        axis.title.x = element_text(size = 16,colour= "#071140FF"), 
        axis.title.y = element_text(size = 16,colour= "#00010DFF"),
        legend.position = "none",
        panel.border = element_rect(colour = "#636363", 
                                    fill=NA,
                                    linewidth=1))

# Save the plot
ggsave(filename = "Cancer Cell line data analysis/figures/WNT2.png",
       plot = wnt2_plot,
       width = 18,
       height = 12,
       units = "in",
       dpi = 300)

# .........................WNT7B.......................................

# Categorize the gene effect scores as Negative or Positive
wnt7b$Gene.effect <- ifelse(wnt7b$Gene.effect.score < 0, 
                                     "negative", "positive")

# Create bar plot
wnt7b_plot <- ggplot(wnt7b, aes(x = Category, 
                                  y = Gene.effect.score,
                                fill = Gene.effect)) +
  geom_bar(stat = "identity", width=0.2) + 
  theme_minimal() +
  labs(title = "Gene Effect Score for WNT7B in Breast Cancer Cell Lines",
       x = "Breast Cancer Cell Lines",
       y = "Gene Effect Score") +
  scale_fill_manual(values = c("negative" = "#9A133DFF", "positive" = "#1A318BFF")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 14, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   colour = "#636363" ),
        axis.text.y = element_text(size = 14,
                                   colour = "#00010DFF" ), 
        plot.title = element_text(size = 20, face = "bold", color = "#071140FF"), 
        axis.title.x = element_text(size = 16,colour= "#071140FF"), 
        axis.title.y = element_text(size = 16,colour= "#00010DFF"),
        legend.position = "none",
        panel.border = element_rect(colour = "#636363", 
                                    fill=NA,
                                    linewidth=1))

# Save the plot
ggsave(filename = "Cancer Cell line data analysis/figures/WNT7B.png",
       plot = wnt7b_plot,
       width = 18,
       height = 12,
       units = "in",
       dpi = 300)



# .........................WNT11.......................................

# Categorize the gene effect scores as Negative or Positive
wnt11$Gene.effect <- ifelse(wnt11$Gene.effect.score < 0, 
                                     "negative", "positive")

# Create bar plot
wnt11_plot <- ggplot(wnt11, aes(x = Category, 
                                y = Gene.effect.score,
                                fill = Gene.effect)) +
  geom_bar(stat = "identity", width=0.2) + 
  theme_minimal() +
  labs(title = "Gene Effect Score for WNT11 in Breast Cancer Cell Lines",
       x = "Breast Cancer Cell Lines",
       y = "Gene Effect Score") +
  scale_fill_manual(values = c("negative" = "#9A133DFF", "positive" = "#1A318BFF")) +
  theme(axis.text.x = element_text(angle = 90, 
                                   size = 14, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   colour = "#636363" ),
        axis.text.y = element_text(size = 14,
                                   colour = "#00010DFF" ), 
        plot.title = element_text(size = 20, face = "bold", color = "#071140FF"), 
        axis.title.x = element_text(size = 16,colour= "#071140FF"), 
        axis.title.y = element_text(size = 16,colour= "#00010DFF"),
        legend.position = "none",
        panel.border = element_rect(colour = "#636363", 
                                    fill=NA,
                                    linewidth=1))

# Save the plot
ggsave(filename = "Cancer Cell line data analysis/figures/WNT11.png",
       plot = wnt11_plot,
       width = 18,
       height = 12,
       units = "in",
       dpi = 300)

# ------------------------------------------------------------------------------
# Combine plots
figure_ccl <- plot_grid(wnt2_plot, wnt7b_plot, wnt11_plot,
                        align = "v",
                        ncol = 1,label_fontfamily = "Open Sans", 
                        label_fontface = "plain",
                        labels = "AUTO",
                        label_size = 40,
                        label_x = -0.01,
                        label_y = 1.02)
# Export the plot
save_plot(filename = "Cancer Cell line data analysis/figures/CCL.png",
          plot = figure_ccl,
          ncol = 1,
          base_height = 18,
          base_width = 22,
          bg = "white")
