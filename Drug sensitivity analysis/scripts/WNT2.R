# Drug Sensitivity Analysis
# Author: Muntasim Fuad

# load required packages 

library(tidyverse)
library(ggpubr)
library(readxl)
library(writexl)

# ...................WNT2............................

# load drug sensitivity data (ctrp)
ctrp_wnt2 <- read_xlsx("Drug sensitivity analysis/data/WNT2/WNT2_CTRP.xlsx")


# Filter data
wnt2_ctrp_neg <- ctrp_wnt2 |> filter(fdr < 0.05) |> 
  filter(cor < -0.12) |> arrange(cor) |> head(20)

# Export top 20 CTRP
write_xlsx(wnt2_ctrp_neg,
           path = "Drug sensitivity analysis/outputs/CTRP/Top 20 CTRP drug sensitivity & WNT2 Expression.xlsx")

# Create the lollipop chart
ctrp_plot2 <- ggplot(wnt2_ctrp_neg, 
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
  labs(title = "Correlation between CTRP drug sensitivity and WNT2 expression",
       x = "Correlation",
       y = "") +
  theme_bw() + theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT2/CTRP.png",
       plot = ctrp_plot2,
       width = 6, height = 7,
       dpi = 300)



# load drug sensitivity data (gdsc)
gdsc_wnt2 <- read_xlsx("Drug sensitivity analysis/data/WNT2/WNT2_GDSC.xlsx")


# Filter data
wnt2_gdsc_neg <- gdsc_wnt2 |> filter(fdr < 0.05) |> 
  filter(cor < -0.08) |> arrange(cor) |> head(10)

wnt2_gdsc_pos <- gdsc_wnt2 |> filter(fdr < 0.1) |> 
  filter(cor > 0.05) |> arrange(cor) |> head(10)

wnt2_gdsc <- bind_rows(wnt2_gdsc_pos, wnt2_gdsc_neg)

# Export top 20 CTRP
write_xlsx(wnt2_gdsc,path = "Drug sensitivity analysis/outputs/GDSC/Top 20 GDSC drug sensitivity & WNT2 Expression.xlsx")

# Create the lollipop chart
gdsc_plot2 <- ggplot(wnt2_gdsc, 
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
  labs(title = "Correlation between GDSC drug sensitivity and WNT2 expression",
       x = "Correlation",
       y = "") +
  theme_bw() + theme(plot.title = element_text(size = 10))

ggsave(filename = "Drug sensitivity analysis/figures/WNT2/GDSC.png",
       plot = gdsc_plot2,
       width = 6, height = 7,
       dpi = 300)
