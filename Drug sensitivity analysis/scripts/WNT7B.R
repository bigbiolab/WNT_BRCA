# Drug Sensitivity Analysis
# Author: Muntasim Fuad

# load required packages 

library(tidyverse)
library(ggpubr)
library(readxl)
library(writexl)

# ...................WNT7B............................

# load drug sensitivity data (ctrp)
ctrp_wnt7b <- read_xlsx("Drug sensitivity analysis/data/WNT7B/WNT7B_CTRP.xlsx")


# Filter data
wnt7b_ctrp_pos <- ctrp_wnt7b |> filter(fdr < 0.05) |> 
  filter(cor > 0.25) |> arrange(cor) |> head(10)

wnt7b_ctrp_neg <- ctrp_wnt7b |> filter(fdr < 0.05) |> 
  filter(cor < -0.12) |> arrange(cor) |> head(10)


wnt7b_ctrp <- bind_rows(wnt7b_ctrp_pos, wnt7b_ctrp_neg)

# Export top 20 CTRP
write_xlsx(wnt7b_ctrp,path = "Drug sensitivity analysis/outputs/CTRP/Top 20 CTRP drug sensitivity & WNT7B Expression.xlsx")

# Create the lollipop chart
ctrp_plot7b <- ggplot(wnt7b_ctrp, 
                    aes(x = cor, y = reorder(drug, cor), 
                        color = cor)) +
  geom_segment(aes(x = 0, xend = cor, 
                   y = reorder(drug, cor), yend = drug), 
               color = "grey") +
  geom_point(aes(size = -log10(fdr)), 
             shape = 20, fill = "black") +
  geom_vline(xintercept = 0, 
             linetype = "solid", 
             color = "grey")+
  scale_color_gradient2(low = "#2b8cbe",
                        mid = "white",
                        high = "#e34a33", midpoint = 0) +
  labs(title = "Correlation between CTRP drug sensitivity and WNT7B expression",
       x = "Correlation",
       y = "") +
  theme_bw() + theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT7B//CTRP.png",
       plot = ctrp_plot7b,
       width = 6, height = 7,
       dpi = 300)



# load drug sensitivity data (gdsc)

gdsc_wnt7b <- read_xlsx("Drug sensitivity analysis/data/WNT7B/WNT7B_GDSC.xlsx")


# Filter data

wnt7b_gdsc_pos <- gdsc_wnt7b |> filter(fdr < 0.1) |> 
  filter(cor > 0.19) |> arrange(cor) |> head(10)

wnt7b_gdsc_neg <- gdsc_wnt7b |> filter(fdr < 0.1) |> 
  filter(cor < -0.2) |> arrange(cor) |> head(10)


wnt7b_gdsc <- bind_rows(wnt7b_gdsc_pos, wnt7b_gdsc_neg)

# Export top 20 CTRP

write_xlsx(wnt7b_gdsc,path = "Drug sensitivity analysis/outputs/GDSC/Top 20 GDSC drug sensitivity & WNT7B Expression.xlsx")

# Create the lollipop chart

gdsc_plot7b <- ggplot(wnt7b_gdsc, 
                    aes(x = cor, y = reorder(drug, cor), 
                        color = cor)) +
  geom_segment(aes(x = 0, xend = cor, 
                   y = reorder(drug, cor), yend = drug), 
               color = "grey") +
  geom_point(aes(size = -log10(fdr)), 
             shape = 20, fill = "black") +
  geom_vline(xintercept = 0, 
             linetype = "solid", 
             color = "grey")+
  scale_color_gradient2(low = "#2b8cbe",
                        mid = "white",
                        high = "#e34a33", midpoint = 0) +
  labs(title = "Correlation between GDSC drug sensitivity and WNT7B expression",
       x = "Correlation",
       y = "") +
  theme_bw() + theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT7B/GDSC.png",
       plot = gdsc_plot7b,
       width = 6, height = 7,
       dpi = 300)
