# GO Enrichment Analysis & Visualization: MF (Molecular Function)  
# Author: Muntasim Fuad

# Load required packages 
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(hrbrthemes)
library(org.Hs.eg.db)
library(openxlsx)

# Load gene data
genes <- read.xlsx("PPIN construction and enrichment analysis/data/genes_string_interaction.xlsx")

# Convert gene data from gene symbol to Entrez ID
genes_enterz <- bitr(genes$Gene.Symbol, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = "org.Hs.eg.db")

genelist <- genes_enterz$ENTREZID

# Perform Gene Ontology (GO) enrichment analysis for Molecular Function (MF)
GO_MF <- enrichGO(
  genelist,
  OrgDb = "org.Hs.eg.db",
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE,
  pool = TRUE
)

# Visualization of GO enrichment result 

# Generate plotting object
mf_bar <- barplot(
  GO_MF,
  x = "Count",
  color = "p.adjust",
  showCategory = 15,
  order = TRUE)

# Customize the plot 
mf_pal <- c("#d95f0e", "#fec44f")

mf_vis <- mf_bar + scale_fill_gradientn(colors = mf_pal)+
  theme_ipsum_rc(base_size = 16, 
                 grid_col = "grey",
                 grid = "XY", 
                 plot_title_margin = 10,
                 axis_title_size = 20,
                 axis_title_just = "c") +
  theme(legend.position="right",
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_text(face = "bold",
                                    size = 24, 
                                    hjust = 0.5,
                                    vjust = 1),
        plot.title = element_text(size = 22, face = "bold")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
  theme(axis.text.y = element_text(margin = margin(r = 5),
                                   size = 20),
        panel.border = element_rect(colour = "#636363", 
                                    fill=NA,
                                    linewidth=1)) +
  ggtitle("Molecular Function (MF)")

# Save the plot
ggsave("PPIN construction and enrichment analysis/figures/MF.png", 
       plot = mf_vis, 
       width = 16, 
       height = 14, 
       dpi = 300,
       units = "in")


# Export Results as table
mf_data <- mf_bar$data
write.csv(mf_data,"PPIN construction and enrichment analysis/outputs/GO_MF.csv")

