
# ======================================================================================
#	1) Import of libraries to manage data frames and definition of arguments
# ======================================================================================

# Import of libraries
library(readr)
library(dplyr)

# Definition of arguments (main_dir refers to the cell type)
args <- commandArgs(trailingOnly = TRUE)
main_path <- args[1]
sample_type <- args[2]

# Definition of the final path.
path <- paste0(main_path, "/", sample_type, "/")



# =============================================================
#       2) Retrieval of the column data from every sample
# =============================================================

sample_type_1_f <- read.table(paste0(path, "subfiles/ColData_<cell>_<sample_type_1>_females.txt"), header = TRUE)
sample_type_2_f <- read.table(paste0(path, "subfiles/ColData_<cell>_<sample_type_2>_females.txt"), header = TRUE)

sample_type_1_m <- read.table(paste0(path, "subfiles/ColData_<cell>_<sample_type_1>_males.txt"), header = TRUE)
sample_type_2_m <- read.table(paste0(path, "subfiles/ColData_<cell>_<sample_type_2>_males.txt"), header = TRUE)



# =========================================================================================
#       3) Combination of the names of the genotype, sample type and sex (identifier)
# =========================================================================================


# 3.1) Identifiers for females.

if(setequal(sample_type_1_f[["genotype"]], sample_type_2_f[["genotype"]])){
	
	unique(sample_type_1_f[["genotype"]])

	genotype_females <- unique(sample_type_1_f[["genotype"]])
	sample_type <- c("<sample_type_1>", "<sample_type_2>")
	sex <- "f"

	comb <- expand.grid(genotype_females, sample_type, sex)
	join_f <- paste(comb$Var1, comb$Var2, comb$Var3, sep = "_")

	print(join_f)

}

# 3.2) Identifiers for females.

if(setequal(sample_type_1_m[["genotype"]], sample_type_2_m[["genotype"]])){

        unique(sample_type_1_m[["genotype"]])

        genotype_males <- unique(sample_type_1_m[["genotype"]])
        sample_type <- c("<sample_type_1>", "<sample_type_2>")
        sex <- "m"

        comb <- expand.grid(genotype_males, sample_type, sex)
        join_m <- paste(comb$Var1, comb$Var2, comb$Var3, sep = "_")

        print(join_m)

}


# =============================================================================================
#       4) Generation of the final output with the combinations of names or identifiers
# =============================================================================================

join_full <- data.frame(c(join_f, join_m))

write_tsv(join_full, paste0(path, "subfiles/names_dir.tsv"))



