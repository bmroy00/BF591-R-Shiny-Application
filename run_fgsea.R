# Extract gene names and corresponding statistics (logFC, p-values, etc.)
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(testthat)
#library(fgsea)

DESeq2_results <- read.csv("data/Diff_exp.csv")
rownames(DESeq2_results) <- DESeq2_results$Column1
deseq2_res_tibble <- as_tibble(DESeq2_results)
deseq2_res_tibble <- deseq2_res_tibble[,-1]

label_res <- function(deseq2_res, padj_threshold) {
  # Convert to tibble
  deseq2_res_tibble <- as_tibble(deseq2_res)
  
  # Create gene column using row names and make it the first column
  deseq2_res_tibble <- mutate(deseq2_res_tibble, genes = rownames(deseq2_res), .before=1)
  
  # Create volc_plot_status col according to log2FoldChange and padj values
  labeled_res <- mutate(
    deseq2_res_tibble, 
    volc_plot_status = case_when(log2FoldChange > 0 & padj < padj_threshold ~ "UP",
                                 log2FoldChange < 0 & padj < padj_threshold  ~ "DOWN",
                                 padj > padj_threshold ~ "NS"
    ), .after=1)
  
  return(labeled_res)
}

labeled_results <- label_res(deseq2_res_tibble, 0.10)
labeled_results <- mutate(labeled_results, genes = symbol)
labeled_results <- labeled_results[,-3]


genes <- DESeq2_results$genes
statistic <- DESeq2_results$log2FoldChange
names(statistic) <- genes

# Create a named vector of statistics (logFC, p-values, etc.)
gene_list <- list(DE = statistic)


if (FALSE) {
# Load gene sets data
pathways <- gmtPathways("m2.cp.v2023.1.Mm.symbols.gmt")

# Run FGSEA
fgsea_results <- fgsea(pathways = pathways, 
             stats = gene_list,
             minSize = 15,
             maxSize = 500)
}