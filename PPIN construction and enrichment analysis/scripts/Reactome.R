# Reactome Enrichment Analysis & Visualization
# Author: Muntasim Fuad

# Load required packages 
library(ReactomePA)
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


# Perform Reactome enrichment analysis
Reactome <- enrichPathway(
  genelist,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = NULL,
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2
)

# Visualization of Reactome Enrichment  result 

# Generate plotting object
reactome_bar <- barplot(
  Reactome,
  x = "Count",
  color = "p.adjust",
  showCategory = 15,
  order = TRUE)

# Customize the plot 
rcpal<- c("#9ca1c3", "#d5ccff","#8856a7")

reactome_vis <- reactome_bar + scale_fill_gradientn(colors = rcpal)+
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
  scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
  theme(axis.text.y = element_text(margin = margin(r = 5),
                                   size = 20),
        panel.border = element_rect(colour = "#636363", 
                                    fill=NA,
                                    linewidth=1)) +
  ggtitle("Reactome")


  
# Save the plot
ggsave("PPIN construction and enrichment analysis/figures/Reactome.png", 
       plot = reactome_vis, 
       width = 18, 
       height = 12, 
       dpi = 300,
       units = "in")


# Export Results as table
reactome_data <- reactome_bar$data
write.csv(reactome_data,"PPIN construction and enrichment analysis/outputs/Reactome.csv")
