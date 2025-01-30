# Differential Gene Expression Analysis: DESeq2
# Author: Md. Jubayer Hossain & Muntasim Fuad

# Load required packages 
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(hgu133a.db)
library(hgu133plus2.db)
library(conflicted)

# Function to download and preprocess data
get_gse_data <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE, 
                destdir = "Validation of prognostic WNTs using GEO/data/")
  if (length(gse) > 1) {
    idx <- grep("GPL", attr(gse, "names"))
    gse <- gse[[idx[1]]]
  } else {
    gse <- gse[[1]]
  }
  counts <- exprs(gse)
  colData <- pData(gse)
  
  # Ensure colData contains a 'condition' column
  if (!"condition" %in% colnames(colData)) {
    colData$condition <- factor(rep(c("control", "treatment"), length.out = nrow(colData)))
  }
  
  # Check if counts are non-integer and round if necessary
  if (!all(counts == round(counts))) {
    counts <- round(counts)
  }
  
  list(counts = counts, colData = colData)
}

# List of GEO dataset IDs
gse_ids <- c("GSE15852", "GSE42568")

# Download and preprocess datasets
datasets <- lapply(gse_ids, get_gse_data)
class(datasets)


# Function to analyze each dataset
analyze_dataset <- function(dataset, id) {
  counts <- dataset$counts
  colData <- dataset$colData
  
  # Check if all counts are zeros
  if (all(counts == 0)) {
    stop("All counts are zero for dataset ", id)
  }
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res$gene <- rownames(res)
  res$study <- id  # Only keep GEO ID
  
  return(as.data.frame(res))
}


results_list <- Map(analyze_dataset, datasets, gse_ids)

# Combine results
combined_results <- bind_rows(results_list)

# Retrive results based on accesion ids
GSE15852 <- combined_results |> filter(study == "GSE15852")
GSE42568 <- combined_results |> filter(study == "GSE42568")

# ------------------------------------------------------------------------------
# Map Probe IDs to gene symbols
GSE15852_symbols <- mapIds(hgu133a.db,
                           keys = GSE15852$gene,
                           column = "SYMBOL",
                           keytype = "PROBEID",
                           multiVals = "first")

GSE42568_symbols <- mapIds(hgu133plus2.db,
                           keys = GSE42568$gene,
                           column = "SYMBOL",
                           keytype = "PROBEID",
                           multiVals = "first")
# ------------------------------------------------------------------------------
# Add gene symbols to the data frame
GSE15852$gene_symbol <- GSE15852_symbols
GSE42568$gene_symbol <- GSE42568_symbols 

# Combine both data frames 
degenes <- bind_rows(GSE15852, GSE42568)

# Rearrange Columns
degenes <-  degenes |> select(gene_symbol,
                             gene,
                             baseMean,
                             log2FoldChange,
                             lfcSE,
                             stat,
                             pvalue,
                             padj,
                             study) |> drop_na(gene_symbol)

# Count the number of genes for each study
gene_counts <- degenes |> 
  group_by(study) |> 
  summarise(number_of_genes = n())

print(gene_counts)

# Export final results
write.csv(degenes,"Validation of prognostic WNTs using GEO/outputs/degenes.csv", 
          row.names = FALSE)