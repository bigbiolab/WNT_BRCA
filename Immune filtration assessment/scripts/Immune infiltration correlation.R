# 1.2.4 Pan-cancer gene expression and immune infiltration correlation
# Author: Shekhar Saha 

# Load required packages
library(TCGAplot)
library(tidyverse)

# Input your gene
gene <- "WNT7B"


# define color gradients
lowcol <- "#3182bd"

highcol <- "#f03b20"


# 1.2.4.1 Pan-cancer gene expression and Immune cell ratio correlation
immucell <-gene_immucell_heatmap(gene,method="pearson",
                                 lowcol = lowcol,
                                 highcol = highcol,
                                 cluster_row=T,
                                 cluster_col=T,legend=T)

ggsave(filename = paste0("Immune filtration assessment/figures/Immune cell ratio correlation/", gene, ".png"), 
       plot = immucell, 
       width = 10,height = 6, 
       dpi = 300)


# 1.2.4.2 Pan-cancer gene expression and Immune score correlation


immunescore <-gene_immunescore_heatmap(gene,method="pearson",
                                       lowcol = lowcol,
                                       highcol = highcol,
                                       cluster_row=T, cluster_col=T,
                                       legend=T)

ggsave(filename = paste0("Immune filtration assessment/figures/Immune score correlation/", gene, ".png"), 
       plot = immunescore, 
       width = 8.5,height = 3, 
       dpi = 300)
