# Cancer Cell Line Analysis
# Author: Muntasim Fuad

# Install required packages
devtools::install_github("gaospecial/ggVennDiagram")

# Load required package
library(ggVennDiagram)
library(tidyverse)
library(grid)

# Load cell line data
WNT2 <- read.csv("Cell line data analysis/data/WNT2.csv")
WNT7B <- read.csv("Cell line data analysis/data/WNT7B.csv")
WNT11 <- read.csv("Cell line data analysis/data/WNT11.csv")

# Filter cell lines negative gene effect score
WNT2 <- WNT2 |> filter(Gene.effect.score < 0) 
WNT7B <- WNT7B |> filter(Gene.effect.score < 0)
WNT11 <- WNT11 |> filter(Gene.effect.score < 0)

# Select cell lines negative gene effect score
WNT2 <- WNT2$Category
WNT7B <- WNT7B$Category
WNT11 <- WNT11$Category

# Create a list of the cell lines
cell_lines_sets <- list(WNT2, WNT7B, WNT11)

# Generate the Venn diagram
ccl_venn <- ggVennDiagram(cell_lines_sets,
                                category.names = c("WNT2", "WNT7B", "WNT11"),
                                set_color = c("#2b8cbe", "#756bb1", "#3182bd"),
                                label_color = "black",
                                label_alpha = 0,
                                label = "count",
                                label_size = 5,
                                set_size = 8) +
  theme(legend.position = "right")+
  scale_fill_gradientn(colors = c("#e7e1ef", "#7fcdbb", "#a1d99b")) +
  labs(title = "Breast Cancer Cell Lines",
       subtitle = "Venn Diagram Visualizing Common Breast Cancer Cell Lines Across Three Genesets",
       fill = "cell lines") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20, 
                                  face = "bold", 
                                  color = "#636363"), 
        plot.subtitle = element_text(size = 15, 
                                     face = "italic", 
                                     color = "#999999"))


# Retrieve intersection of three gene sets
intersect <- process_region_data(Venn(cell_lines_sets))

intersect_ccl <- intersect |> 
  filter(name == "Set_1/Set_2/Set_3") |> 
  pull(item) |>  as.data.frame()

colnames(intersect_ccl) <- "Breast Cancer Cell Lines"

# Convert intersect_genes to a string for annotation
ccl_text <- intersect_ccl$`Breast Cancer Cell Lines`

# Combine the text with newline characters for multiline display
ccl_text_combined <- paste(ccl_text, collapse = "\n")

# Create the text grob (text graphic object)
ccl_grob <- textGrob(ccl_text_combined, x = 0.84,
                     y = 0.04,  # Adjust this value as needed
                     hjust = 0, vjust = 0,
                     gp = gpar(fontsize = 16))

# Create the rectangle grob (for the frame) with adjusted position
rect_grob <- rectGrob(x = unit(1, "npc") - unit(1.5, "lines"),
                      y = unit(0, "npc") + unit(1.5, "lines"),
                      width = grobWidth(ccl_grob) + unit(1.5, "lines"),
                      height = grobHeight(ccl_grob) + unit(1.5, "lines"),
                      just = c("right", "bottom"),
                      gp = gpar(col = "#023743FF", fill = "#E7E9E4FF"))

combined_grob <- grobTree(rect_grob, ccl_grob)

# Save the plot
ggsave(filename = "Cell line data analysis/figures/Venn diagram.png",
       plot = ccl_venn +
         annotation_custom(combined_grob) +
         geom_point(aes(x = 2, 
                        y = -2), 
                    size = 3, 
                    color = "#023743FF") +
         geom_segment(aes(x = 2, 
                          y = -2, 
                          xend = 6.337, 
                          yend = -5.362),
                      colour = "#023743FF",
                      linewidth = 1),
       width = 12,
       height = 12,
       dpi = 300, 
       units = "in")
