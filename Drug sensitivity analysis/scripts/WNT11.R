# Drug Sensitivity Analysis
# Author: Muntasim Fuad

# load required packages 

library(tidyverse)
library(ggpubr)
library(readxl)
library(writexl)

# ...................WNT11............................

# load drug sensitivity data (ctrp)

ctrp_wnt11 <- read_xlsx("Drug sensitivity analysis/data/WNT11/WNT11_CTRP.xlsx")


# Filter data

wnt11_ctrp_pos <- ctrp_wnt11 |> filter(fdr < 0.1) |> 
  filter(cor > 0.001) |> arrange(cor) |> head(10)

wnt11_ctrp_neg <- ctrp_wnt11 |> filter(fdr < 0.1) |> 
  filter(cor < -0.12) |> arrange(cor) |> head(10)


wnt11_ctrp <- bind_rows(wnt11_ctrp_pos, wnt11_ctrp_neg)

# Export top 20 CTRP

write_xlsx(wnt11_ctrp,path = "Drug sensitivity analysis/outputs/CTRP/Top 20 CTRP drug sensitivity & WNT11 Expression.xlsx")

# Create the lollipop chart

ctrp_plot11 <- ggplot(wnt11_ctrp, 
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
  labs(title = "Correlation between CTRP drug sensitivity and WNT11 expression",
       x = "Correlation",
       y = "") +
  theme_bw() + theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT11//CTRP.png",
       plot = ctrp_plot11,
       width = 6, height = 7,
       dpi = 300)



# load drug sensitivity data (gdsc)

gdsc_wnt11 <- read_xlsx("Drug sensitivity analysis/data/WNT11/WNT11_GDSC.xlsx")


# Filter data

wnt11_gdsc_pos <- gdsc_wnt11 |> filter(fdr < 0.05) |> 
  filter(cor > 0.11) |> arrange(cor) |> head(10)

wnt11_gdsc_neg <- gdsc_wnt11 |> filter(fdr < 0.1) |> 
  filter(cor < -0.04) |> arrange(cor) |> head(10)


wnt11_gdsc <- bind_rows(wnt11_gdsc_pos, wnt11_gdsc_neg)

# Export top 20 CTRP

write_xlsx(wnt11_gdsc,path = "Drug sensitivity analysis/outputs/GDSC/Top 20 GDSC drug sensitivity & WNT11 Expression.xlsx")

# Create the lollipop chart

gdsc_plot11 <- ggplot(wnt11_gdsc, 
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
  labs(title = "Correlation between GDSC drug sensitivity and WNT11 expression",
       x = "Correlation",
       y = "") +
  theme_bw() +theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT11/GDSC.png",
       plot = gdsc_plot11,
       width = 6, height = 7,
       dpi = 300)
