# KEGG Enrichment Analysis & Visualization
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


# Perform KEGG enrichment analysis
KEGG <- enrichKEGG(
  genelist,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)

# Visualization of KEGG enrichment result 

# Generate plotting object
kegg_bar <- barplot(
  KEGG,
  x = "Count",
  color = "p.adjust",
  showCategory = 15,
  order = TRUE)

# Customize the plot 
kcpal<- c("#3182bd", "#9ecae1")

kegg_vis <- kegg_bar + scale_fill_gradientn(colors = kcpal)+
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
                                    size = 22, 
                                    hjust = 0.5,
                                    vjust = 1),
        plot.title = element_text(size = 22, face = "bold")) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 60)) +
          theme(axis.text.y = element_text(margin = margin(r = 5),
                                           size = 20),
                panel.border = element_rect(colour = "#636363", 
                                            fill=NA,
                                            linewidth=1)) +
  ggtitle("KEGG")
 
# Save the plot        
ggsave("PPIN construction and enrichment analysis/figures/KEGG.png", 
       plot = kegg_vis, 
       width = 16, 
       height = 10, 
       dpi = 300,
       units = "in")


# Export Results as table
kegg_data <- kegg_bar$data
write.csv(kegg_data,"PPIN construction and enrichment analysis/outputs/KEGG.csv")
