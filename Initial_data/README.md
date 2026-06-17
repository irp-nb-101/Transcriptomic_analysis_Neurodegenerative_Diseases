
# Transcriptomics methods to study neurodegenerative diseases (still in development)


## About the project

Differential expression analysis (DEA) with bulk and cell-specific data. The comparison is carried out between controls and different genotypes. The main steps are:

- Filtering of the data frames.
- Establishment of a PCA to select the design.
- Performing of DEA with DESeq2.



## Environment setup

The tools required to carry out the analysis need to be downloaded through Bioconductor in R:

```r
BiocManager::install(c("DESeq2"), force = TRUE)
```

```r
BiocManager::install(c("EnhancedVolcano"), force = TRUE)
```

```r
BiocManager::install("IHW")
```



## Scripts execution: retrieval of the results

The performance of (DEA), as well as the corresponding representations with Volcano Plot, is carried out in the following script: RScript_DESeq_AD.R.

Regarding the fields that need to be modified, here they are listed and explained:

- path <- "" (line 21): the directory where the matrix counts and coldata are located.
- males_csv <- paste0(path, "") (line 27): the matrix counts of males.
- females_csv <- paste0(path, "") (line 28): the matrix counts of females.
- males_info <- paste0(path, "") (line 45): coldata of males.
- females_info <- paste0(path, "") (line 46): coldata of females.

- colnames_sel <- ... grepl("", ... (line 87): search for the string in the columns that match the rows of the coldata.
- alter_genes <- c("") (line 122): manual addition of genes that might add noise to the analysis.
- design = ~ (lines 153 and 204): design of the experiment (one applied before the PCA and another after the PCA).
- pcaData <- ...intgroup = c("", ...; ...color = ...; ...distance_matrix ~ ... (lines 156-174 and 208-220): condition (different than "Control") used for the PCA.

- path_new <- "" (line 198): path where all the generated data frames with the log2FC and padj-values will be stored.
- path_viz <- "" (line 277): path where all the generated Volcano plots will be stored.


Final execution (if the script and the initial files are stored in Linux):

```r
Rscript RScript_DESeq_AD.R
```


## Results (final output)

The main result is based on a Volcano plot for every condition tested in the experiment for both males and females (separately).

