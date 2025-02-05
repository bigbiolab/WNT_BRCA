# Single-Cell Analysis
# Author: Muntasim Fuad

# Load required package
library(ggVennDiagram)
library(tidyverse)
library(grid)

# Load single-cell data
WNT2 <- read.csv("Single cell analysis/data/TISCH_WNT2_heatmap.csv")
WNT7B <- read.csv("Single cell analysis/data/TISCH_WNT7B_heatmap.csv")
WNT11 <- read.csv("Single cell analysis/data/TISCH_WNT11_heatmap.csv")

# Filter cell types with expression values > 0
WNT2 <- WNT2 |> filter(Values > 0) 
WNT7B <- WNT7B |> filter(Values > 0)
WNT11 <- WNT11 |> filter(Values > 0)

# Select cell lines negative gene effect score
WNT2 <- WNT2$Category
WNT7B <- WNT7B$Category
WNT11 <- WNT11$Category

# Create a list of the cell lines
single_cell_sets <- list(WNT2, WNT7B, WNT11)

# Generate the Venn diagram
ccl_venn <- ggVennDiagram(single_cell_sets,
                                category.names = c("WNT2", "WNT7B", "WNT11"),
                                set_color = c("#756bb1", "#c51b8a", "#43a2ca"),
                                label_color = "black",
                                label_alpha = 0,
                                label = "count",
                                label_size = 5,
                                set_size = 8) +
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c("#fee0d2", "#fc9272", "#e34a33")) +
  labs(title = "Single-Cell Analysis",
       subtitle = "Venn Diagram Visualizing Common Cell Types Across Three Genesets",
       fill = "cell lines") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20, 
                                  face = "bold", 
                                  color = "#636363"), 
        plot.subtitle = element_text(size = 15, 
                                     face = "italic", 
                                     color = "#999999"))


# Retrieve intersection of three gene sets
intersect <- process_region_data(Venn(single_cell_sets))

intersect_sc <- intersect |> 
  filter(name == "Set_1/Set_2/Set_3") |> 
  pull(item) |>  as.data.frame()

colnames(intersect_sc) <- "Cell Types"

# Convert intersect_genes to a string for annotation
sc_text <- intersect_sc$`Cell Types`

# Combine the text with newline characters for multiline display
sc_text_combined <- paste(sc_text, collapse = "\n")

# Create the text grob (text graphic object)
sc_grob <- textGrob(sc_text_combined, x = 0.82,
                     y = 0.04,  # Adjust this value as needed
                     hjust = 0, vjust = 0,
                     gp = gpar(fontsize = 16))

# Create the rectangle grob (for the frame) with adjusted position
rect_grob <- rectGrob(x = unit(1, "npc") - unit(1.5, "lines"),
                      y = unit(0, "npc") + unit(1.5, "lines"),
                      width = grobWidth(sc_grob) + unit(1.5, "lines"),
                      height = grobHeight(sc_grob) + unit(1.5, "lines"),
                      just = c("right", "bottom"),
                      gp = gpar(col = "#023743FF", fill = "#E7E9E4FF"))

combined_grob <- grobTree(rect_grob, sc_grob)

# Save the plot
ggsave(filename = "Single cell analysis/figures/Venn diagram.png",
       plot = ccl_venn +
         annotation_custom(combined_grob) +
         geom_point(aes(x = 2, 
                        y = -2), 
                    size = 3, 
                    color = "#023743FF") +
         geom_segment(aes(x = 2, 
                          y = -2, 
                          xend = 6.02, 
                          yend = -5.77),
                      colour = "#023743FF",
                      linewidth = 1),
       width = 12,
       height = 12,
       dpi = 300, 
       units = "in",
       bg = "white")
