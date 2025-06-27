# Differential Gene Expression Analysis
# Author: Muntasim Fuad

# 00. Load required packages
library(tidyverse)
library(DESeq2)
library(biomaRt)
library(EnhancedVolcano)
library(openxlsx)

conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(base::as.factor)
conflicted::conflicts_prefer(dplyr::rename)

# 01. Load raw counts data
# GEO Accession: GSE233242
GSE233242 <- read_tsv("data/GSE233242_raw_counts_GRCh38.p13_NCBI.tsv")
glimpse(GSE233242)

# 02. Load metadata
meta_data <- read_csv("data/GSE233242_metadata.csv")
glimpse(meta_data)

# 03. Create a matrix & add gene ids as row names
count_data <- GSE233242 |> 
  column_to_rownames("GeneID") |>
  as.matrix()

# 04. Match metadata with count data
meta_data <- meta_data |> 
  filter(Sample %in% colnames(count_data)) |> 
  arrange(match(Sample, colnames(count_data)))

# 05. Prepare Sample information
colData <- data.frame(condition = as.factor(meta_data$Condition), 
                      row.names = colnames(count_data))


# 06. Create DESeq2 data set object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

# 07. Filter out lowly expressed genes
# Keep genes with at least 10 counts in at least half the samples
keep_genes <- rowSums(count_data >= 10) >= (ncol(count_data) / 2)
count_data <- count_data[keep_genes, ]

# 08. Run DESeq2 analysis
dds <- DESeq(dds)

# 09. Extract results
res <- results(dds)
res$Gene_ID <- rownames(res)
res$Study <- "GSE233242"

# 10. Set up biomaRt for annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

# 11. Fetch gene annotations (NCBI Gene ID → Symbol/Description)
annotations <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "description"),
  filters = "entrezgene_id", 
  values = res$Gene_ID,
  mart = mart
)

# 12. Merge annotations directly with the combined results
annotated_results <- merge(as.data.frame(res), 
                           annotations, by.x = "Gene_ID", 
                           by.y = "entrezgene_id", 
                           all.x = TRUE)

# 13. Rename columns for clarity
annotated_results <- annotated_results |> 
  rename(Gene_Symbol = hgnc_symbol, 
         Gene_Description = description)

# 14. Save results to a CSV file
write.csv(annotated_results, "outputs/GSE233242_DE_results.csv")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 01. Load raw counts data
# GEO Accession: GSE227679
GSE227679 <- read_tsv("data/GSE227679_raw_counts_GRCh38.p13_NCBI.tsv")
glimpse(GSE227679)

# 02. Load metadata
meta_data <- read_csv("data/GSE227679_metadata.csv")
glimpse(meta_data)

# 03. Create a matrix & add gene ids as row names
count_data <- GSE227679 |> 
  column_to_rownames("GeneID") |>
  as.matrix()

# 04. Match metadata with count data
meta_data <- meta_data |> 
  filter(Sample %in% colnames(count_data)) |> 
  arrange(match(Sample, colnames(count_data)))

# 05. Prepare Sample information
colData <- data.frame(condition = as.factor(meta_data$Condition), 
                      row.names = colnames(count_data))


# 06. Create DESeq2 data set object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

# 07. Filter out lowly expressed genes
# Keep genes with at least 10 counts in at least half the samples
keep_genes <- rowSums(count_data >= 10) >= (ncol(count_data) / 2)
count_data <- count_data[keep_genes, ]

# 08. Run DESeq2 analysis
dds <- DESeq(dds)

# 09. Extract results
res <- results(dds)
res$Gene_ID <- rownames(res)
res$Study <- "GSE227679"

# 10. Set up biomaRt for annotation
mart <- useMart("ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")

# 11. Fetch gene annotations (NCBI Gene ID → Symbol/Description)
annotations <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "description"),
  filters = "entrezgene_id", 
  values = res$Gene_ID,
  mart = mart
)

# 12. Merge annotations directly with the combined results
annotated_results <- merge(as.data.frame(res), 
                           annotations, by.x = "Gene_ID", 
                           by.y = "entrezgene_id", 
                           all.x = TRUE)

# 13. Rename columns for clarity
annotated_results <- annotated_results |> 
  rename(Gene_Symbol = hgnc_symbol, 
         Gene_Description = description)

# 14. Save results to a CSV file
write.csv(annotated_results, "outputs/GSE227679_DE_results.csv")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Volcano plot for visualization of differentially expressed genes
# Author: Muntasim Fuad

# Load data
GSE227679_DE_results <- read_csv("outputs/GSE227679_DE_results.csv")
GSE233242_DE_results <- read_csv("outputs/GSE233242_DE_results.csv")

# Combine DE results from both studies
degenes <- bind_rows(GSE227679_DE_results, GSE233242_DE_results) |>
  select(Gene_ID, Gene_Symbol,Gene_Description, log2FoldChange, padj, Study) |> 
  drop_na(Gene_Symbol)

# Generate volcano plot
volcano <- EnhancedVolcano(degenes,
                           lab = as.character(degenes$Gene_Symbol),
                           x = 'log2FoldChange',
                           y = 'padj',
                           xlim = c(-10, 10),
                           pointSize = 2.0,
                           labSize = 6.0,
                           col=c('#636363', '#9ecae1', '#a1d99b', '#e6550d'),
                           colAlpha = 1,
                           legendLabels=c('NS','Log2FC','adj.p-value',
                                          'adj.p-Value & Log2FC'),
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
ggsave(filename = ("figures/Volcano Plot.png"), 
       plot = volcano, 
       width = 9, 
       height = 8, 
       dpi = 600)

# Export combined results as .xlsx file
write.xlsx(degenes, "outputs/merged_DE_results.xlsx")

# ------------------------------------------------------------------------------
# Count the number of genes for each study
gene_counts <- degenes |> 
  group_by(Study) |> 
  summarise(number_of_genes = n())

print(gene_counts)

# Count differentially expressed genes
GSE227679_DEGs <- degenes |> 
  filter(padj < 0.05 & Study == "GSE227679") |> 
  summarise(GSE227679 = n())

GSE227679_DEGs

GSE233242_DEGs <- degenes |> 
  filter(padj < 0.05 & Study == "GSE233242") |> 
  summarise(GSE233242 = n())

GSE233242_DEGs

Total_DEGs <- degenes |> filter(padj < 0.05)|> 
  summarise(number_of_genes = n())

Total_DEGs

# Calculate sample number for each study 
meta_data |> count(Condition)
