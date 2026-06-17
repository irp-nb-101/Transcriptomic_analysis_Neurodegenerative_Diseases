
#========================================================================================================
#   1) Loading of DESeq2 and additional packages to perform the DEA + Addition of the main arguments
#========================================================================================================


# Loading of all the useful packages.

library(dplyr)
library(readr)
library(ggplot2)
library(utils)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)
library(IHW)
library(vegan)


# Definition of the main arguments of the execution

args <- commandArgs(trailingOnly = TRUE)
cell_pop <- args[1]
path <- args[2]
sample_type <- args[3]
pattern <- args[4]


# Manual modification of the genes that might add noise to the analysis.

alter_genes <- c("")



#======================================================================================================
#   2) Retrieval and necessary modifications of the count matrices and their information data
#======================================================================================================


# Storage of the data frames containing the count matrix.
males_csv <- paste0(path, "matrix_counts/matrix_males_", cell_pop, ".csv")
females_csv <- paste0(path, "matrix_counts/matrix_females_", cell_pop, ".csv")
males_df_c <- read.csv(males_csv)
males_df_def_c <- males_df_c %>%
  select(-width) # (*) The column "width" might be named differently in each case.
females_df_c <- read.csv(females_csv)
females_df_def_c <- females_df_c %>%
  select(-width) # (*) The column "width" might be named differently in each case.

# Joining of the females count matrix to the males one.
full_df_values <- females_df_c %>%
  full_join(males_df_c, by = "X")

full_df_values_def <- females_df_def_c %>%
  full_join(males_df_def_c, by = "X")


# Storage of the data frames containing the information of each sample (cohort, sex, genotype and training).
males_info <- paste0(path, "col_data/ColData_", cell_pop, "_", sample_type, "_males.txt")
females_info <- paste0(path, "col_data/ColData_", cell_pop, "_", sample_type, "_females.txt")
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

colnames_sel <- colnames(full_df_values)[grepl(pattern, colnames(full_df_values))]

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


com_1 <- alter_genes[!(alter_genes %in% common_noisy_genes)]
com_1

common_total <- c(common_noisy_genes, com_1)

females_final_filt_df <- females_df_values_filt_def %>%
  filter(!(rownames(females_df_values_filt_def) %in% common_total))

males_final_filt_df <- males_df_values_filt_def %>%
  filter(!(rownames(males_df_values_filt_def) %in% common_total))




#======================================================================================
#   4) Automation of the Principal Component Analysis (PCA) - Covariables selection
#======================================================================================


# 4.1) Automation in females.

# Definition of variables (with the change of position to establish the order).

vars <- colnames(females_info_df_def)[colnames(females_info_df_def) %in% c("training", "IP.cohort")]

if("IP.cohort" %in% vars){
  tmp1 <- vars[1]
  tmp2 <- vars[2]
  vars[1] <- tmp2
  vars[2] <- tmp1
  
  # Creation of the first DESeq2 design with all the variables.
  
  des_f <- DESeqDataSetFromMatrix(countData = females_final_filt_df, colData = females_info_df_def, design = ~ IP.cohort + training + genotype)
} else {
  des_f <- DESeqDataSetFromMatrix(countData = females_final_filt_df, colData = females_info_df_def, design = ~ training + genotype)
}

# Normalisation of the DESeq2 object.

norm_des_f <- DESeq(des_f)
vsd <- vst(norm_des_f, blind = TRUE)


# Performing of a PCA with each one of the covariables, including a PERMANOVA test to determine significant differences.

noise_vars <- c()
for(v in vars){
  
  pcaData <- plotPCA(vsd, intgroup = c(v, "genotype"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca <- ggplot(pcaData, aes(PC1, PC2, color = .data[[v]], shape = genotype)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed()
  
  ggsave(paste0(path, "PCA/pca_", v, "_", cell_pop, "_", sample_type, "_f.png"), 
         plot = pca, 
         width = 10,  
         height = 10,  
         units = "in", 
         dpi = 300)
  
  distance_matrix <- dist(t(assay(vsd)))
  
  design_vector <- as.formula(
    paste("distance_matrix ~", v, "+ genotype")
  )
  
  adonis_res <- adonis2(design_vector, data = as.data.frame(colData(vsd)))
  
  if(adonis_res[["Pr(>F)"]][1] < 0.05){
    noise_vars <- c(noise_vars, v)
  }
  
}

covariables <- paste(noise_vars, collapse = " + ")
new_design_f <- as.formula(paste0("~ ", covariables, " + genotype"))
environment(new_design_f) <- globalenv()



# 4.2) Automation in males.

# Definition of variables (with the change of position to establish the order).

vars <- colnames(males_info_df_def)[colnames(males_info_df_def) %in% c("training", "IP.cohort")]

if("IP.cohort" %in% vars){
  tmp1 <- vars[1]
  tmp2 <- vars[2]
  vars[1] <- tmp2
  vars[2] <- tmp1
  
  # Creation of the first DESeq2 design with all the variables.
  
  des_m <- DESeqDataSetFromMatrix(countData = males_final_filt_df, colData = males_info_df_def, design = ~ IP.cohort + training + genotype)
} else {
  des_m <- DESeqDataSetFromMatrix(countData = males_final_filt_df, colData = males_info_df_def, design = ~ training + genotype)
}


# Normalisation of the DESeq2 object.

norm_des_m <- DESeq(des_m)
vsd <- vst(norm_des_m, blind = TRUE)


# Performing of a PCA with each one of the covariables, including a PERMANOVA test to determine significant differences.

noise_vars <- c()
for(v in vars){
  
  pcaData <- plotPCA(vsd, intgroup = c(v, "genotype"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca <- ggplot(pcaData, aes(PC1, PC2, color = .data[[v]], shape = genotype)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed()
  
  ggsave(paste0(path, "PCA/pca_", v, "_", cell_pop, "_", sample_type, "_m.png"), 
         plot = pca, 
         width = 10,  
         height = 10,  
         units = "in", 
         dpi = 300)
  
  distance_matrix <- dist(t(assay(vsd)))
  
  design_vector <- as.formula(
    paste("distance_matrix ~", v, "+ genotype")
  )
  
  adonis_res <- adonis2(design_vector, data = as.data.frame(colData(vsd)))
  
  if(adonis_res[["Pr(>F)"]][1] < 0.05){
    noise_vars <- c(noise_vars, v)
  }
  
}

covariables <- paste(noise_vars, collapse = " + ")
new_design_m <- as.formula(paste0("~ ", covariables, " + genotype"))
environment(new_design_m) <- globalenv()




#========================================================================================================
#   5) Performing of the Differential Expression Analysis for both sexes with the selected covariables
#========================================================================================================


# 5.1) Differential Expression Analysis in females.

des_f <- DESeqDataSetFromMatrix(countData = females_final_filt_df, colData = females_info_df_def, design = new_design_f)
norm_des_f <- DESeq(des_f)

levels_gen <- levels(norm_des_f$genotype)[levels(norm_des_f$genotype) != "Control"]
control <- "Control"


# Retrieval of the results from all the possible combinations of controls and samples with the specified genotypes.

results_vs_control <- lapply(levels_gen, function(p) {
  results(norm_des_f, contrast=c("genotype", p, "Control"), filterFun = ihw, alpha = 0.05)
})

names(results_vs_control) <- paste0(levels_gen, "_vs_", control)
path_new <- paste0(path, "results/DESeq2_res_", cell_pop, "_", sample_type, "/")


# Storage of the results.

for(i in 1:length(results_vs_control)){
  res_final_f <- as.data.frame(results_vs_control[[i]])
  res_final_sorted_f <- res_final_f %>%
    arrange(padj)
  genotype <- gsub("_vs_Control", "", names(results_vs_control)[i])
  write.table(res_final_sorted_f, paste0(path_new, "dea_Ctrl-vs-", genotype, "-", sample_type, "_f.tsv"), sep = "\t")
}



# 5.2) Differential Expression Analysis in males.

des_m <- DESeqDataSetFromMatrix(countData = males_final_filt_df, colData = males_info_df_def, design = new_design_m )
norm_des_m <- DESeq(des_m)

levels_gen <- levels(norm_des_m$genotype)[levels(norm_des_m$genotype) != "Control"]
control <- "Control"


# Retrieval of the results from all the possible combinations of controls and samples with the specified genotypes.

results_vs_control_m <- lapply(levels_gen, function(p) {
  results(norm_des_m, contrast=c("genotype", p, "Control"), filterFun = ihw, alpha = 0.05)
})

names(results_vs_control_m) <- paste0(levels_gen, "_vs_", control)


# Storage of the results.

for(i in 1:length(results_vs_control_m)){
  res_final_m <- as.data.frame(results_vs_control_m[[i]])
  res_final_sorted_m <- res_final_m %>%
    arrange(padj)
  genotype <- gsub("_vs_Control", "", names(results_vs_control_m)[i])
  write.table(res_final_sorted_m, paste0(path_new, "dea_Ctrl-vs-", genotype, "-", sample_type, "_m.tsv"), sep = "\t")
}


