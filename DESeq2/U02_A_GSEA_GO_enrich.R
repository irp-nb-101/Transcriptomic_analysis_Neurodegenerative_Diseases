
#=================================================================================================
#   1) Loading of the necessary packages to perform the GSEA + Addition of the main arguments
#=================================================================================================


# Loading of all the useful packages.

library(utils)
library(DESeq2)
library(apeglm)
library(EnhancedVolcano)
library(RColorBrewer)
library(ashr)
library(BiocManager)
library(sva)
library(org.Mm.eg.db)
library(DOSE)
library(AnnotationHub)
library(ensembldb)
library(clusterProfiler)
library(ggplot2)
library(gtable)
library(pheatmap)
library(readr)
library(dplyr)
library(stringr)
library(gprofiler2)


# Definition of the main arguments of the execution.

args <- commandArgs(trailingOnly = TRUE)
path_new_files <- args[1]
dir <- args[2]


# Listing of the files from the DESeq2 results directory.

list_files <- list.files(paste0(path_new_files, "results/DESeq2_", dir, "/"))



#==============================================================
#   2) Perfoming of the Gene Set Enrichment Analysis (GSEA)
#==============================================================

for(file in list_files){
  if(str_detect(file, "dea_")){
    
    # Creation of the suffix of the file for further storage.
    name_suffix_file <- str_split((str_split(file, "dea_")[[1]][2]), "\\.tsv")[[1]][1]
    
    # GSEA execution through the function gost() from gprofiler2.
    matrix_not_sorted <- read.table(paste0(path_new_files, "results/DESeq2_", dir, "/", file), header = TRUE)
    matrix_sorted <- matrix_not_sorted %>%
      arrange(stat)
    gostres_def <- gost(query = rownames(matrix_sorted), sources = c("GO:BP", "KEGG", "REAC", "WP"), organism = "mmusculus", ordered_query = TRUE)
    gostres_def_df <- as.data.frame(gostres_def[[1]])
    
    
    if(length(rownames(gostres_def_df)) > 0){
    
    gostres_def_df_sort <- gostres_def_df[1:13] %>%
      arrange(term_name)
    write.table(gostres_def_df_sort, paste0(path_new_files, "results/GSEA_", dir, "/gsea_", name_suffix_file, ".tsv"), sep = "\t")
  
    res <- gostres_def_df %>%
      head(30)
   
    res_reg <- gostres_def_df %>%
      arrange(substr(term_name, 1, 1) < "M", term_name) %>%
      head(30)
    
  	res$ratio <- res$intersection_size / res$query_size
  	res_reg$ratio <- res_reg$intersection_size / res_reg$query_size
  	
  	
  	# Enrichment plot representation of the global results.
  
  	pl <- ggplot(res, aes(x = ratio, y = reorder(term_name, ratio), size = intersection_size, color = p_value)) +
    		geom_point() +
    		scale_color_gradient(low = "red", high = "blue") +
    		theme_minimal() +
  	    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  	    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  	  theme(
  	    axis.text.y = element_text(size = 7.5, lineheight = 0.6), 
  	    panel.grid.major.y = element_blank() 
  	  )
  	
  	
  	# Enrichment plot representation of the results based on regulation terms.
  	
  	pl_reg <- ggplot(res_reg, aes(x = ratio, y = reorder(term_name, ratio), size = intersection_size, color = p_value)) +
  	  geom_point() +
  	  scale_color_gradient(low = "orange", high = "purple") +
  	  theme_minimal() +
  	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  	  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  	  theme(
  	    axis.text.y = element_text(size = 7.5, lineheight = 0.6), 
  	    panel.grid.major.y = element_blank()
  	  )
  
  	ggsave(paste0(path_new_files, "results/GSEA_", dir, "/gsea_", name_suffix_file, "_plot.jpg"), pl)
  	ggsave(paste0(path_new_files, "results/GSEA_", dir, "/gsea_", name_suffix_file, "_plot_reg.jpg"), pl_reg)
  	
	}
  }
}



#========================================================
#   3) Perfoming of the Gene Ontology (GO) enrichment
#========================================================


for(file in list_files){
  if(str_detect(file, "dea_")){
    
    # Creation of the suffix of the file for further storage.
    name_suffix_file <- str_split((str_split(file, "dea_")[[1]][2]), "\\.tsv")[[1]][1]
    
    # Import of the DESeq2 results matrix and definition of the significant genes (global, downregulated and upregulated).
    
    matrix_not_sorted <- read.table(paste0(path_new_files, "results/DESeq2_", dir, "/", file), header = TRUE)
    all_genes_input <- as.character(rownames(matrix_not_sorted))
    res <- as.data.frame(matrix_not_sorted)

    signif_res_input <- res[
      res$padj < 0.05 & !is.na(res$padj),
    ]

    signif_genes_input <- as.character(rownames(signif_res_input))

    signif_res_input_down <- res[
      res$padj < 0.05 &
        !is.na(res$padj) &
        res$log2FoldChange < 0,
    ]

    signif_res_input_up <- res[
      res$padj < 0.05 &
        !is.na(res$padj) &
        res$log2FoldChange > 0,
    ]

    signif_genes_input_up <- as.character(rownames(signif_res_input_up))
    signif_genes_input_down <- as.character(rownames(signif_res_input_down))
  

    # GO analysis with the global significant genes, including the corresponding enrichment plots.
    
    print(length(signif_genes_input))
    if(length(signif_genes_input) >= 10) {

      ego_BP_input <- enrichGO(gene = signif_genes_input, universe = all_genes_input,
                                        keyType = "SYMBOL",
                                        OrgDb = org.Mm.eg.db,
                                        ont = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        readable = TRUE)
      
      cluster_summary_BP <- as.data.frame(ego_BP_input)

      cluster_readable <- setReadable(ego_BP_input, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
      cluster_readable_df <- as.data.frame(cluster_readable)

      write.table(cluster_readable_df, paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, ".tsv"), sep = "\t")

      
      cluster_readable_df_num <- cluster_readable_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        head(30)
      
      cluster_readable_df_reg <- cluster_readable_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        arrange(substr(Description, 1, 1) < "M", Description) %>%
        head(30)
      
      
      pl <- ggplot(cluster_readable_df_num, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(x = "Gene ratio") +
        labs(y = "Pathway name") +
        scale_color_gradient(low = "red", high = "blue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), # Ajusta el espacio entre líneas rotas
          panel.grid.major.y = element_blank() # Limpia el gráfico para facilitar la lectura
        )
      
      
      pl_reg <- ggplot(cluster_readable_df_reg, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(y = "Pathway name") +
        labs(x = "Gene ratio") +
        scale_color_gradient(low = "orange", high = "purple") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), # Ajusta el espacio entre líneas rotas
          panel.grid.major.y = element_blank() # Limpia el gráfico para facilitar la lectura
        )
      

	    ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_plot.jpg"), pl)
	    ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_plot_reg.jpg"), pl_reg)
	    
    }


    
    # GO analysis with the significant upregulated genes, including the corresponding enrichment plots.
    
    print(length(signif_genes_input_up))
    if(length(signif_genes_input_up) >= 10) {

      ego_BP_input_up <- enrichGO(gene = signif_genes_input_up, universe = all_genes_input,
                                        keyType = "SYMBOL",
                                        OrgDb = org.Mm.eg.db,
                                        ont = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        readable = TRUE)
      
      cluster_summary_BP_up <- as.data.frame(ego_BP_input_up)

      cluster_readable_up <- setReadable(ego_BP_input_up, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
      cluster_readable_up_df <- as.data.frame(cluster_readable_up)

      write.table(cluster_readable_up_df, paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_up.tsv"), sep = "\t")
      
      cluster_readable_df_num <- cluster_readable_up_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        head(30)
      
      cluster_readable_df_reg <- cluster_readable_up_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        arrange(substr(Description, 1, 1) < "M", Description) %>%
        head(30)
      
      
      pl <- ggplot(cluster_readable_df_num, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(x = "Gene ratio") +
        labs(y = "Pathway name") +
        scale_color_gradient(low = "red", high = "blue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), # Ajusta el espacio entre líneas rotas
          panel.grid.major.y = element_blank() # Limpia el gráfico para facilitar la lectura
        )
      
      
      pl_reg <- ggplot(cluster_readable_df_reg, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(y = "Pathway name") +
        labs(x = "Gene ratio") +
        scale_color_gradient(low = "orange", high = "purple") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), # Ajusta el espacio entre líneas rotas
          panel.grid.major.y = element_blank() # Limpia el gráfico para facilitar la lectura
        )
      
      ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_up_plot.jpg"), pl)
      ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_up_plot_reg.jpg"), pl_reg)
    }

    
    
    # GO analysis with the significant downregulated genes, including the corresponding enrichment plots.

    print(length(signif_genes_input_down))
    if(length(signif_genes_input_down) >= 10) {
      ego_BP_input_down <- enrichGO(gene = signif_genes_input_down, universe = all_genes_input,
                                        keyType = "SYMBOL",
                                        OrgDb = org.Mm.eg.db,
                                        ont = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.05,
                                        readable = TRUE)
      
      cluster_summary_BP_down <- as.data.frame(ego_BP_input_down)

      cluster_readable_down <- setReadable(ego_BP_input_down, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")
      cluster_readable_down_df <- as.data.frame(cluster_readable_down)

      write.table(cluster_readable_down_df, paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_down.tsv"), sep = "\t")
    
      
      cluster_readable_df_num <- cluster_readable_down_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        head(30)
      
      cluster_readable_df_reg <- cluster_readable_down_df %>%
        mutate(GeneRatio_number = sapply(parse(text = GeneRatio), eval)) %>%
        arrange(substr(Description, 1, 1) < "M", Description) %>%
        head(30)
      
      
      pl <- ggplot(cluster_readable_df_num, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(x = "Gene ratio") +
        labs(y = "Pathway name") +
        scale_color_gradient(low = "red", high = "blue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), 
          panel.grid.major.y = element_blank()
        )
      
      
      pl_reg <- ggplot(cluster_readable_df_reg, aes(x = GeneRatio_number, y = reorder(Description, GeneRatio_number), size = Count, color = qvalue)) +
        geom_point() +
        labs(y = "Pathway name") +
        labs(x = "Gene ratio") +
        scale_color_gradient(low = "orange", high = "purple") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        theme(
          axis.text.y = element_text(size = 7.5, lineheight = 0.6), 
          panel.grid.major.y = element_blank() 
        )
      
      ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_down_plot.jpg"), pl)
      ggsave(paste0(path_new_files, "results/GO_", dir, "/go_", name_suffix_file, "_down_plot_reg.jpg"), pl_reg)
      
    }
  }
}



