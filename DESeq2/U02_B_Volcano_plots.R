
#========================================================================================================
#   1) Loading of DESeq2 and additional packages to perform the DEA + Addition of the main arguments
#========================================================================================================


# Loading of all the useful packages.

library(dplyr)
library(readr)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(stringr)
library(EnhancedVolcano)                                                                                                                                                                                                                                                                                                


# Definition of the main arguments of the execution.

args <- commandArgs(trailingOnly = TRUE)

path_new_files <- args[1]
dir <- args[2]


# Listing of the files from the DESeq2 results directory.

list_files <- list.files(paste0(path_new_files, "results/DESeq2_", dir, "/"))



#==============================================================
#   2) Volcano plots creation for each DESeq2 results table 
#==============================================================


for(file in list_files){
  if(str_detect(file, "dea_")){
    
    # Creation of the suffix of the file for further storage.
    name_suffix_file <- str_split((str_split(file, "dea_")[[1]][2]), "\\.tsv")[[1]][1]
    
    # Import of DESeq2 results table and corresponding filtering.
    
    df_ex <- read.table(paste0(path_new_files, "results/DESeq2_", dir, "/", file), header = TRUE)
  
    df_ex_sign <- df_ex %>%
      mutate(significant = padj < 0.05 & abs(log2FoldChange) >= 2) %>%
      arrange(padj)

    # Retrieval of the 15 most significant genes for each parameter and labels definition.
    
    ordered_p_value <- df_ex_sign %>%
      arrange(padj) %>%
      head(15)
    
    ordered_log2FC <- df_ex_sign %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(15)
  
    df_ex$label <- ""
    df_ex$label[rownames(df_ex) %in% rownames(ordered_p_value) | rownames(df_ex) %in% rownames(ordered_log2FC)] <- 
      rownames(df_ex)[rownames(df_ex) %in% rownames(ordered_p_value) | rownames(df_ex) %in% rownames(ordered_log2FC)]
    df_ex_sign_labels <- rownames(df_ex_sign)
    
    
    # Creation of the volcano plot with the specified labels.

    EV <- EnhancedVolcano(df_ex,
                  title = "Enhanced Volcano with Airways",
                  lab = df_ex$label,
                  selectLab = rownames(df_ex_sign_labels),
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlab = 'log2FoldChange',
                  ylab = '-log10(FDR)',
                  pCutoff = 0.05,
                  FCcutoff = 2,
                  max.overlaps = Inf,
                  drawConnectors = TRUE)

    
    ggsave(paste0(path_new_files, "results/Volcano_", dir, "/VP_", name_suffix_file, ".png"),
       plot = EV,
       width = 10,
       height = 10,
       units = "in",
       dpi = 300)

  }
}

