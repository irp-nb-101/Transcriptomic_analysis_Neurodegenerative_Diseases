

# ======================================================================================
#       1) Import of libraries to manage data frames and definition of arguments
# ======================================================================================


library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

main_path <- args[1]
cell_pop <- args[2]
dir <- args[3]

print(main_path)
print(cell_pop)
print(dir)

path <- paste0(main_path, "/", cell_pop, "/")
print(path)


# ================================================================================================================
#       2) Import of the column data file and the BAM ID identifiers file according to each directory (dir)
# ================================================================================================================


if(grepl("<sample_type_1>_f", dir)){

info <- paste0(path, "subfiles/ColData_<cell>_<sample_type_1>_females.txt")
info_df <- read.table(info)

bam <- paste0(path, "subfiles/File_Sample_<cell>_females.txt")
bam_df <- read.table(bam)

}


if(grepl("<sample_type_2>_f", dir)){

info <- paste0(path, "subfiles/ColData_<cell>_<sample_type_2>_females.txt")
info_df <- read.table(info)

bam <- paste0(path, "subfiles/File_Sample_<cell>_females.txt")
bam_df <- read.table(bam)

}


if(grepl("<sample_type_1>_m", dir)) {

info <- paste0(path, "subfiles/ColData_<cell>_<sample_type_1>_males.txt")
info_df <- read.table(info)

bam <- paste0(path, "subfiles/File_Sample_<cell>_males.txt")
bam_df <- read.table(bam)

}


if(grepl("<sample_type_2>_m", dir)){

info <- paste0(path, "subfiles/ColData_<cell>_<sample_type_2>_males.txt")
info_df <- read.table(info)

bam <- paste0(path, "subfiles/File_Sample_<cell>_males.txt")
bam_df <- read.table(bam)


}

cols <- c()
for(i in info_df[1,]){
  cols <- c(cols, i)
}
colnames(info_df) <- cols

info_df_def <- info_df %>%
  filter(genotype != "genotype")


cols <- c()
for(i in bam_df[1,]){
  cols <- c(cols, i)
}
colnames(bam_df) <- cols

bam_df_def <- bam_df %>%
  filter(SampleName != "SampleName") # The name of the column is replaced for each case.

colnames(bam_df_def)[3] <- "sample"



# ================================================================================================================
#       3) Filter the number of the short identifiers of each file according to the corresponding genotype
# ================================================================================================================


bam_filt_df_def <- bam_df_def %>%
  semi_join(info_df_def, by = "sample")

bam_filt2_df_def <- bam_filt_df_def %>%
  select(-FileName2)

for(i in 1:length(bam_filt2_df_def[[1]])){
  identifier <- (strsplit(bam_filt2_df_def[[1]][i], "_"))[[1]][2]
  bam_filt2_df_def[[1]][i] <- identifier
}


final_def_df <- info_df_def %>%
  full_join(bam_filt2_df_def, by = "sample")

write_tsv(final_def_df, paste0(path, "subfiles/complete_tables/complete_", dir, ".tsv"))


pattern <- "^(.*)_(.*)_(.*)$"
matches <- regmatches(dir, regexec(pattern, dir))
genotype_select <- matches[[1]][2]


filtered_def_df <- final_def_df %>%
  select(genotype, FileName1) %>%  # The name of the column "FileName1" is replaced for each case.
  filter(genotype == genotype_select) %>%
  select(FileName1)


write_tsv(filtered_def_df, paste0(path, "subfiles/filtered_bam_names/filt_", dir, ".tsv"))




