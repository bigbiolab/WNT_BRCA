# Volcano plot for visualization of differentially expressed genes
# Author: Muntasim Fuad

# Load required pacakges
library(EnhancedVolcano)
library(tidyverse)

# Load data 
degenes <- read.csv("Validation of prognostic WNTs using GEO/outputs/degenes.csv")

# Generate volcano plot
volcano <- EnhancedVolcano(degenes,
                           lab = as.character(degenes$gene_symbol),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           pointSize = 2.0,
                           labSize = 5.0,
                           col=c('#636363', '#9ecae1', '#a1d99b', '#e6550d'),
                           colAlpha = 1,
                           legendLabels=c('NS','Log2FC','p-value',
                                          'p-Value & Log2FC'),
                           legendPosition = 'top',
                           legendLabSize = 16,
                           legendIconSize = 5.0,
                           title = "",  
                           titleLabSize = 16,
                           subtitle = "",
                           subtitleLabSize = 18,
                           pCutoff = 10e-6,
                           FCcutoff = 2,
                           cutoffLineType = "dashed",
                           border = "partial")

# Export volcano plot
ggsave(filename = ("Validation of prognostic WNTs using GEO/figures/vocano.png"), 
       plot = volcano, 
       width = 8.5, 
       height = 8, 
       dpi = 300)