# Gene-Gene Correlation Analysis using TCGAplot
# Author: Muntasim Fuad

# Load required packages
library(TCGAplot)
library(tidyverse)

# Input your 1st gene symbol
gene1 <- "WNT2"

# Input your 2nd gene symbol
gene2 <- "WNT7B"

# Input your 2nd gene symbol
gene3 <- "WNT11"

# Input cancer type
cancer <- "BRCA"

# Generate Gene-gene correlation scatter plot 
png(filename = paste0("figures/Gene-gene correlation scatter plot/", gene1, "_", gene2, ".png"), 
    width = 6000, height = 5000, 
    res = 600)

gene_gene_scatter(cancer,gene1,gene2,density="T") # whether density of gene expression was shown.

dev.off()

# ------------------------------------------------------------------------------
png(filename = paste0("figures/Gene-gene correlation scatter plot/", gene2, "_", gene3, ".png"), 
    width = 6000, height = 5000, 
    res = 600)

gene_gene_scatter(cancer,gene2,gene3,density="T") # whether density of gene expression was shown.

dev.off()
# ------------------------------------------------------------------------------
png(filename = paste0("figures/Gene-gene correlation scatter plot/", gene1, "_", gene3, ".png"), 
    width = 6000, height = 5000, 
    res = 600)

gene_gene_scatter(cancer,gene1,gene3,density="T") # whether density of gene expression was shown.

dev.off()
