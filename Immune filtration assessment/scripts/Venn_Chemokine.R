# Chemokines
# Author: Muntasim Fuad

# Load required package
library(ggVennDiagram)
library(grid)
library(tidyverse)
library(readxl)

# Load data
chemokine <- read_xlsx("Immune filtration assessment/data/Chemokine correlation.xlsx")
  
  
# Define the gene lists
WNT2 <- chemokine$WNT2
WNT7B <- chemokine$WNT7B
WNT11 <- chemokine$WNT11

# Create a list of the gene sets
gene_sets <- list(WNT2, WNT7B, WNT11)

# Generate the Venn diagram
chemo_venn <- ggVennDiagram(gene_sets,
                                category.names = c("WNT2", "WNT7B", "WNT11"),
                                set_color = c("#f03b20", "#3182bd", "#31a354"),
                                label_color = "black",
                                label_alpha = 0,
                                label = "count",
                                label_size = 5,
                                set_size = 8) +
  scale_fill_gradientn(colors = c("#efedf5", "#bcbddc", "#756bb1")) +
  labs(title = "Chemokines",
       subtitle = "Venn Diagram Visualizing Common Chemokines Across Three Genesets",
       fill = "chemokines") +
  theme(legend.position = "right", 
        plot.title = element_text(size = 20, 
                                  face = "bold", 
                                  color = "#636363"), 
        plot.subtitle = element_text(size = 15, 
                                     face = "italic", 
                                     color = "#999999"))

# Retrieve intersection of three gene sets
intersect <- process_region_data(Venn(gene_sets))

intersect_genes <- intersect |> 
  filter(name == "Set_1/Set_2/Set_3") |> 
  pull(item) |>  as.data.frame()

colnames(intersect_genes) <- "chemokines"

# Convert intersect_genes to a string for annotation
gene_text <- paste(intersect_genes$chemokines, collapse = "\n")

# Create the text grob (text graphic object)
gene_grob <- textGrob(gene_text, x = unit(1, "npc") - unit(1, "lines"),
                      y = unit(0, "npc") + unit(1, "lines"),
                      hjust = 1.25, vjust = -0.11,
                      gp = gpar(fontsize = 16))

# Create the rectangle grob (for the frame) with adjusted position
rect_grob <- rectGrob(x = unit(1, "npc") - unit(1.5, "lines"),
                      y = unit(0, "npc") + unit(1.5, "lines"),
                      width = grobWidth(gene_grob) + unit(1.5, "lines"),
                      height = grobHeight(gene_grob) + unit(1.5, "lines"),
                      just = c("right", "bottom"),
                      gp = gpar(col = "#023743FF", fill = "#E7E9E4FF"))

combined_grob <- grobTree(rect_grob, gene_grob)


# Save the plot
ggsave(filename = "Immune filtration assessment/figures/Chemokines.png",
       plot = chemo_venn + 
         annotation_custom(combined_grob) +
         geom_point(aes(x = 2, 
                        y = -2), 
                    size = 3, 
                    color = "#023743FF") +
         geom_segment(aes(x = 2, 
                          y = -2, 
                          xend = 6.8, 
                          yend = -6.58),
                      colour = "#023743FF",
                      linewidth = 1),
       width = 12,
       height = 12,
       dpi = 300, 
       units = "in")

