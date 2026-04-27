

#=======================================================================
#   1) Loading of DESeq2 and additional packages to perform the DEA
#=======================================================================

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(IHW)



#======================================================================================================
#   2) Retrieval and necessary modifications of the count matrices and their information data
#======================================================================================================


# Definition of the path that leads to the files.
path <- ""

# Loading of the dplyr package to manage data frames.
library(dplyr)

# Storage of the data frames containing the count matrix.
males_csv <- paste0(path, "")
females_csv <- paste0(path, "")
males_df_c <- read.csv(males_csv)
males_df_def_c <- males_df_c %>%
  select(-width)
females_df_c <- read.csv(females_csv)
females_df_def_c <- females_df_c %>%
  select(-width)

# Joining of the females count matrix to the males one.
full_df_values <- females_df_c %>%
  full_join(males_df_c, by = "X")

full_df_values_def <- females_df_def_c %>%
  full_join(males_df_def_c, by = "X")


# Storage of the data frames containing the information of each sample (cohort, sex, genotype and training).
males_info <- paste0(path, "")
females_info <- paste0(path, "")
males_info_df <- read.table(males_info)
females_info_df <- read.table(females_info)


# Establish the names of the first row (the real colnames) as the colnames.
cols <- c()
for(i in males_info_df[1,]){
  cols <- c(cols, i)
}
colnames(males_info_df) <- cols

cols <- c()
for(i in females_info_df[1,]){
  cols <- c(cols, i)
}
colnames(females_info_df) <- cols


# Removal of the first row (redundant) and addition of a new column that represents the corresponding sex.
males_info_df_def <- males_info_df %>%
  filter(genotype != "genotype") %>%
  mutate(sex = rep("male", length(males_info_df[,1]) - 1))

females_info_df_def <- females_info_df %>%
  filter(genotype != "genotype") %>%
  mutate(sex = rep("female", length(females_info_df[,1]) - 1))


# Joining of the rows from both data frames (females and males).
full_info_df_def <- bind_rows(females_info_df_def, males_info_df_def)



#=============================================================================================================
#   3) Filtering the common values in both data frames and the genes that lead to analysis artifacts.
#=============================================================================================================


# 3.1) Filtering the column of the values data frame that are only present as rows in the info data frame.

colnames_sel <- colnames(full_df_values)[grepl("", colnames(full_df_values))]

full_df_values_filt <- full_df_values %>%
  select(X, colnames_sel) %>%
  select(X, full_info_df_def[[1]])

rownames(full_df_values_filt) <- full_df_values_filt[,1]
full_df_values_filt_def <- full_df_values_filt %>%
  select(-X)


# 3.2) Filtering the the rows with a total sum of less than 10 counts across the female samples.

females_df_values_filt_def <- full_df_values_filt_def[colnames(full_df_values_filt_def) %in% females_info_df_def[[1]]]

females_df_values_filt_inf <- rowSums(females_df_values_filt_def) < 10
females_df_values_filt_names_df <- females_df_values_filt_def[females_df_values_filt_inf,]
females_df_values_filt_names <- rownames(females_df_values_filt_names_df)


# 3.3) Filtering the rows with a total sum of less than 10 counts across the male samples.

males_df_values_filt_def <- full_df_values_filt_def[colnames(full_df_values_filt_def) %in% males_info_df_def[[1]]]

males_df_values_filt_inf <- rowSums(males_df_values_filt_def) < 10
males_df_values_filt_names_df <- males_df_values_filt_def[males_df_values_filt_inf,]
males_df_values_filt_names <- rownames(males_df_values_filt_names_df)


# 3.4) Filtering both females and males data frames regarding the noisy genes found in both data frames. 

common_noisy_genes_logical <- females_df_values_filt_names %in% males_df_values_filt_names
common_noisy_genes_logical
common_noisy_genes <- females_df_values_filt_names[females_df_values_filt_names %in% males_df_values_filt_names]

alter_genes <- c("")  # Inclusion of additional genes that might add noise to the analysis.

com_1 <- alter_genes[!(alter_genes %in% common_noisy_genes)]
com_1

common_total <- c(common_noisy_genes, com_1)
length(common_noisy_genes)
length(com_1)
length(common_total)

females_final_filt_df <- females_df_values_filt_def %>%
  filter(!(rownames(females_df_values_filt_def) %in% common_total))

males_final_filt_df <- males_df_values_filt_def %>%
  filter(!(rownames(males_df_values_filt_def) %in% common_total))




#================================================================================================================
#   4) Performing the RNA-Seq analyses in both sexes (including the PCA to determine the design in each case)
#================================================================================================================


# 4.1) RNA-Seq analysis in females.


# A + C) Initial design (before PCA) and updated design (after PCA and possible modification of the design).

library(utils)

des_f <- DESeqDataSetFromMatrix(countData = females_final_filt_df, colData = females_info_df_def, design = ~ ) # With all possible variables.
norm_des_f <- DESeq(des_f)

vsd <- vst(norm_des_f, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup = c("", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


# B) PCA to determine the variability between the genotype and an additional variable.

library(ggplot2)

ggplot(pcaData, aes(PC1, PC2, color = , shape = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()

library(vegan)

distance_matrix <- dist(t(assay(vsd)))
adonis_res <- adonis2(distance_matrix ~, data = as.data.frame(colData(vsd)))


# D) Retrieval of the results regarding the expression comparison between the control and each modified genotype.

levels_gen <- levels(norm_des_f$genotype)[levels(norm_des_f$genotype) != "Control"]
control <- "Control"

results_vs_control <- lapply(levels_gen, function(p) {
  results(norm_des_f, contrast=c("genotype", p, "Control"), filterFun = ihw, alpha = 0.05)
})

names(results_vs_control) <- paste0(levels_gen, "_vs_", control)


# E) Storage of results.

path_new <- ""

for(i in 1:length(results_vs_control)){
  res_final_f <- as.data.frame(results_vs_control[[i]])
  res_final_sorted_f <- res_final_f %>%
    arrange(padj)
  write.csv(res_final_sorted_f, paste0(path_new, "dea_", names(results_vs_control)[i], "_f.tsv"), sep="\t")
}



# 4.2) RNA-Seq analysis in males.

des_m <- DESeqDataSetFromMatrix(countData = males_final_filt_df, colData = males_info_df_def, design = ~ )
norm_des_m <- DESeq(des_m)


vsd <- vst(norm_des_m, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup = c("", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = , shape = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed()

distance_matrix <- dist(t(assay(vsd)))
adonis_res <- adonis2(distance_matrix ~, data = as.data.frame(colData(vsd)))



levels_gen <- levels(norm_des_m$genotype)[levels(norm_des_m$genotype) != "Control"]

control <- "Control"


results_vs_control_m <- lapply(levels_gen, function(p) {
  results(norm_des_m, contrast=c("genotype", p, "Control"), filterFun = ihw, alpha = 0.05)
})

names(results_vs_control_m) <- paste0(levels_gen, "_vs_", control)


for(i in 1:length(results_vs_control_m)){
  res_final_m <- as.data.frame(results_vs_control_m[[i]])
  res_final_sorted_m <- res_final_m %>%
    arrange(padj)
  write.csv(res_final_sorted_m, paste0(path_new, "dea_", names(results_vs_control_m)[i], "_m.tsv"), sep="\t")
}




#===========================================================
#   5) Representation of the results in a Volcano Plot
#===========================================================


# Iterative construct that stores the Volcano Plot representation for each file created in step 4.

files_path <- list.files(path_new)

for(file in files_path){
  
  print(file)
  
  df_ex <- read.csv(paste0(path, file))
  
  gene_list <- df_ex %>%
    select(2) %>%
    tibble::rownames_to_column("ID")
  
  df_ex_sign <- df_ex %>%
    mutate(significant = padj < 0.05 & abs(log2FoldChange) >= 2) %>%
    arrange(padj)
  
  top_genes <- df_ex_sign$X[c(1:6)]
  
  EV <- EnhancedVolcano(df_ex_sign,
                        title = "",
                        lab = df_ex$X,
                        x = 'log2FoldChange',
                        y = 'padj')
  
  path_viz <- ""
  
  ggsave(paste0(path_viz, "VP_", file,".png"), 
         plot = EV, 
         width = 10,  
         height = 10,  
         units = "in", 
         dpi = 300)
  
}




