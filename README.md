
# Transcriptomics methods to study neurodegenerative diseases (still in development)


## About the project

Differential expression analysis (DEA) with bulk and cell-specific data. The comparison is carried out between controls and different genotypes. The main steps are:

- Filtering of the data frames.
- Establishment of a PCA to select the design.
- Performing of DEA with DESeq2.
- Creation of the GSEA and GO enrichment plots
- Quantitative visualisation with Volcano plots.



## Environment setup (with some of the examples of the packages listed in the scripts)

The tools required to carry out the analysis need to be downloaded through Bioconductor in R:

```r
BiocManager::install("DESeq2")
```

```r
BiocManager::install("IHW")
```

```r
BiocManager::install("gprofiler2")
```

```r
BiocManager::install("EnhancedVolcano")
```



## Scripts execution: retrieval of the results

The performance of (DEA), as well as the corresponding representations with Volcano Plot, is carried out in the following script: RScript_DESeq_AD.R.

Regarding the fields that need to be modified in the subscript U01_A_RScript_DESeq_AD.R, here they are listed and explained:
- alter_genes <- c("") (line 31): manual addition of genes that might add noise to the analysis.

As for the fields that need to be modified in the main script S01_DESeq2_dea.slurm, here they are also listed:
- arg 1 --> cell_pop (cell population)
- arg 2 --> $PWD/ (current directory)
- arg 3 --> sample_type
- arg 4 --> pattern

- Rscript U01_A_RScript_DESeq_AD.R # {1} {2} {3} {4} --> Modification site inside of the script


### DESeq2 DEA script execution through SLURM

```bash
sbatch S01_DESeq2_dea.slurm
```

### Creation of the Volcano and enrichment plots.

Afterwards, the script S02_Graph_dea.slurm consists of the execution through SLURM job arrays of the R scripts that create the graphical results for the differential expression analysis.

One modification must be carried out: the number of job arrays (minus 1) in the job directive --array=0-X, where X indicates the total number of elements of an array (n) minus 1 (X = n - 1).

```bash
sbatch S02_Graph_dea.slurm
```


## Results (final output)

The main results obtained for every condition tested in the experiment for both males and females (separately) are based on enrichment plots with the pathways associated to the differentially expressed genes and a Volcano plot that shows the level of significance of these genes regarding the False Discovery Rate (FDR) value and log2Fchange (Log2FC) or expression coefficient between conditions. 




# Introduction to the alternative splicing analysis approach in neurodegenerative diseases (still in development).


## About the project

Multivariate Analysis of Transcript Splicing (rMATS) that considers all the different types of splicing events. The currently available steps are:

- Formatting of the BAM files according to the selected conditions for each comparison.
- Execution of rMATS.


## Environment setup (rMATS installment in a Conda environment)

```bash
conda create -n rmats python=3.10
```

```bash
conda config --env --add channels conda-forge
```

```bash
conda config --env --add channels bioconda
```

```bash
conda install rmats
```


## Location of the scripts used for this part of the analysis.

The scripts that are currently available in the Alternative_Splicing folder are based on the main rMATS execution script S03_rmats.slurm (in the rMATS subfolder) and the previous formatting scripts (in the subscripts subfolder).


### Examples of execution with previously indicated variables (through the use of --export).

```bash
sbatch --export=CELL=<cell_name> S02_format_files.slurm
```

```bash
sbatch --export=SP=<single/paired>,L=<read_length> S03_rmats.slurm
```


