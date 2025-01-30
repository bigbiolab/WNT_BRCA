# Pan-cancer gene expression and immune-related genes correlation
# Author: Muntasim Fuad

# Load required packages

library(TCGAplot)
library(tidyverse)


# Input your gene symbol
gene <- "WNT11"

# define color gradients
lowcol <- "#3182bd"

highcol <- "#f03b20"

# 1.2.3.1 Pan-cancer gene expression and ICGs correlation
checkpoint <- gene_checkpoint_heatmap(gene,
                                      method="pearson",
                                      lowcol= lowcol,
                                      highcol= highcol,
                                      cluster_row=T,
                                      cluster_col=T,
                                      legend=T)

ggsave(filename = paste0("Immune filtration assessment/figures/ICGs (immune checkpoint genes) correalation/", gene, ".png"), 
       plot = checkpoint,
       width = 8, height = 3, 
       dpi = 300)



# 1.2.3.2 Pan-cancer gene expression and Chemokine correlation 

chemokine <- gene_chemokine_heatmap(gene,
                                    method="pearson",
                                    lowcol= lowcol,
                                    highcol= highcol,
                                    cluster_row=T,
                                    cluster_col=T,
                                    legend=T)
  
ggsave(filename = paste0("Immune filtration assessment/figures/Chemokine correlation/", gene, ".png"),plot = chemokine, 
       width = 8, height = 8,
       dpi = 300)


# 1.2.3.3 Pan-cancer gene expression and Chemokine receptor correlation

receptor <- gene_receptor_heatmap(gene,
                                  method="pearson",
                                  lowcol= lowcol,
                                  highcol= highcol,
                                  cluster_row=T,cluster_col=T,
                                  legend=T)


ggsave(filename = paste0("Immune filtration assessment/figures/Chemokine receptor correlation/", gene, ".png"), 
       plot = receptor, 
       width = 8, 
       height = 5,
       dpi = 300)


# 1.2.3.4 Pan-cancer gene expression and Immune stimulator correlation

immustimulator <-gene_immustimulator_heatmap(gene,
                                             method="pearson",
                                             lowcol= lowcol,
                                             highcol= highcol,
                                             cluster_row=T,
                                             cluster_col=T,
                                             legend=T)


ggsave(filename = paste0("Immune filtration assessment/figures/Immune stimulator correlation/", gene, ".png"), 
       plot = immustimulator, 
       width = 8,height = 9, 
       dpi = 300)


# 1.2.3.5 Pan-cancer gene expression and Immune inhibitor correlation

immuinhibitor <- gene_immuinhibitor_heatmap(
  gene,method="pearson",
  lowcol= lowcol, highcol= highcol,
  cluster_row=T, cluster_col=T,
  legend=T
)

ggsave(filename = paste0("Immune filtration assessment/figures/Immune inhibitor correlation/", gene, ".png"), 
       plot = immuinhibitor, 
       width = 8,height = 6, 
       dpi = 600)
