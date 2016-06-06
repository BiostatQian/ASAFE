## ------------------------------------------------------------------------
# Clear workspace and load ASAFE
rm(list=ls())
library(ASAFE)

# Rows: Marker ID
# Cols: Marker ID, 2 consecutive columns for each individual's 2 chromosomes 
# (i.e. ADM1, ADM1, ADM2, ADM2, ..., ADM250, ADM250)
dim(adm_ancestries) # 56,003 x 501

# Rows: Marker ID
# Cols: Marker ID, 1 column for each individual
# (i.e. ADM1, ADM2, ..., ADM250)
dim(adm_genotypes) # 56003 x 251


## ------------------------------------------------------------------------

# Making the rsID column row names

row.names(adm_ancestries) <- adm_ancestries[,1]
row.names(adm_genotypes) <- adm_genotypes[,1]

adm_ancestries <- adm_ancestries[,-1]
adm_genotypes <- adm_genotypes[,-1]


## ------------------------------------------------------------------------

# apply() works on each row (SNP) of the genotypes matrix, and returns a list,
# with each element of the list corresponding to a SNP.
#
# alleles_list is a list of lists. 
# Elements of the outer list correspond to SNPs.
# Elements of the inner list correspond to 250 
# individual's alleles with no delimiter "/" separating alleles. 
alleles_list = apply(X = adm_genotypes, MARGIN = 1, FUN = strsplit, split = "/")

# Creates a matrix: Number of alleles (ADM1, ADM1, ..., ADM250, ADM250) x (SNPs)
alleles_unlisted = sapply(alleles_list, unlist)

# Change elements of the matrix to numeric:
# Number of alleles (ADM1, ADM1, ..., ADM250, ADM250) x (SNPs).
alleles = apply(X = alleles_unlisted, MARGIN = 2, as.numeric)


## ------------------------------------------------------------------------
# Make Table 1 of the paper
#
# - Rows: 
# 1) Anc 0 Mean Error, 2) Anc 0 SD Error
# 3) Anc 1 Mean Error, 4) Anc 1 SD Error
# 5) Anc 2 Mean Error, 6) Anc 2 SD Error
#
# - Cols: True Allele 1 frequency bins
# 1) (0, 0.2], 2) (0.2, 0.4], 3) (0.4, 0.6], 
# 4) (0.6, 0.8], 5) (0.8, 1.0]

results = matrix(nrow = 6, ncol = 5)

# Get mean, sd errors for ancestry 0, and all 5 true allele freuqency bins
results[1:2,] = get_mean_sd_error_anc(ancestry = 0, 
                                      estimates = adm_estimates, 
                                      truth = ancestral_freqs)

# Get mean, sd errors for ancestry 1, and all 5 true allele freuqency bins
results[3:4,] = get_mean_sd_error_anc(ancestry = 1, 
                                      estimates = adm_estimates, 
                                      truth = ancestral_freqs)

# Get mean, sd errors for ancestry 2, and all 5 true allele freuqency bins
results[5:6,] = get_mean_sd_error_anc(ancestry = 2, 
                                      estimates = adm_estimates, 
                                      truth = ancestral_freqs)

results

## ------------------------------------------------------------------------

n_ind <- 250
n_markers <- 5000

# Load in functions needed to use the ASAFE EM algorithm
library(ASAFE)

# Source in all functions needed to assess error
# when [p0, p1, p2] takes different values, and
# when there is error in local ancestry calls.

source("../R/draw_allele_given_anc.R")
source("../R/get_true_freqs_1snp.R")
source("../R/get_errors_1_scenario.R")
source("../R/get_errors_summary_stats_1_scenario.R")
source("../R/get_scenario_errors.R")
source("../R/sample_ancestry.R")
source("../R/change_ancestry.R")
source("../R/get_results_error.R")

# Ancestry-specific allele frequencies.
# Rows: Scenarios.
# Cols: P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
# Each row is a [P(Allele 1 | Anc 0, P(Allele 1 | Anc 1, P(Allele 1 | Anc 2)]
# combination that we will consider.

anc_spec_freqs <- matrix(c(0, 0.1, 0.4, 0.48, 0.1, 0.4,
		         0.5, 0.5, 0.5, 0.5, 0.6, 0.5,
			 0.9, 0.8, 0.6, 0.52, 0.6, 0.5),
			 nrow = 6, ncol = 3)

#### Draw ancestries

set.seed(1)

# Ancestries matrix.
# Rows: Individuals' alleles. Cols: Markers.

ancestries_matrix <- matrix(NA, nrow = 2 * n_ind, ncol = n_markers)

# Randomly generate ancestries
# for each individual's allele at every marker.
# Ancestries are 0, 1, or 2.

for(snp in 1:n_markers){

    ancestries_matrix[,snp] <- replicate(n = 2 * n_ind, 
                                         sample(x = 0:2, size = 1))
    
}

# Throw out the draws (i.e. columns i.e. markers)
# where there are not all three ancestries (0,1, or 2).

col_all_3_ancs <- apply(  X = ancestries_matrix,
	                      MARGIN = 2,
	                      FUN = function(col_vector){

	                          ind_all_3 <- (length(unique(col_vector)) == 3)
	                          return(ind_all_3)

                          })

# This fixed ancestries matrix will be used for multiple
# ancestral allele frequency scenarios.

ancestries_matrix_original <- ancestries_matrix[, col_all_3_ancs]

# Error rate = 0. out_0.0 should be similar to Supplementary Table 1.

error_rate <- 0

out_0.0 <- get_results_error(error_rate = error_rate,
                             anc_spec_freqs = anc_spec_freqs,
                             ancestries_matrix_true = ancestries_matrix_original)

out_0.0

# write.table(x = out_0.0, file = paste("error_summary_stats_error_", error_rate,
# 			 "_jobID_", jobID,
# 			 ".txt", sep = ""),
#                          quote = FALSE,
#                          col.names = TRUE, row.names = TRUE, append = FALSE)

# Error rate = 0.07. out_0.07 should be similar to Supplementary Table 2.

error_rate <- 0.07

out_0.07 <- get_results_error(error_rate = error_rate,
                         anc_spec_freqs = anc_spec_freqs,
                         ancestries_matrix_true = ancestries_matrix_original)

out_0.07

# write.table(x = out_0.07, file = paste("error_summary_stats_error_", error_rate,
#                         "_jobID_", jobID,
#                         ".txt", sep = ""),
#                         quote = FALSE,
#                         col.names = TRUE, row.names = TRUE, append = FALSE)


## ------------------------------------------------------------------------
# Clear workspace and load ASAFE
rm(list=ls())
library(ASAFE)

# adm_ancestries_test is a matrix with
# Rows: Markers
# Columns: Marker ID, individuals' chromosomes' ancestries
# (e.g. ADM1, ADM1, ADM2, ADM2, and etc.)

# adm_genotypes_test is a matrix with
# Rows: Markers
# Columns: Marker ID, individuals' genotypes (a1/a2)
# (e.g. ADM1, ADM2, ADM3, and etc.)

adm_ancestries_test <- head(adm_ancestries)
adm_genotypes_test <- head(adm_genotypes)

# Making the rsID column row names

row.names(adm_ancestries_test) <- adm_ancestries_test[,1]
row.names(adm_genotypes_test) <- adm_genotypes_test[,1]

adm_ancestries_test <- adm_ancestries_test[,-1]
adm_genotypes_test <- adm_genotypes_test[,-1]

# alleles_list is a list of lists.
# Outer list elements correspond to SNPs.
# Inner list elements correspond to 250 people's alleles 
# with no delimiter separating alleles.
alleles_list <- apply(X = adm_genotypes_test, MARGIN = 1, 
                      FUN = strsplit, split = "/")

# Creates a matrix: 
# Alleles for chromosomes (ADM1, ADM1, ..., ADM250, ADM250) x (SNPs)
alleles_unlisted <- sapply(alleles_list, unlist)

# Change elements of the matrix to numeric
alleles <- apply(X = alleles_unlisted, MARGIN = 2, as.numeric)

# Apply the EM algorithm to each SNP to obtain
# ancestry-specific allele frequency estimates for all SNPs in
# matrices alleles and adm_ancestries_test.
#
# Columns correspond to markers.
# Rows correspond to ancestries 0, 1, and then 2.
# Entries in rows 2 through 4
# give P(Allele 1 | Ancestry a), a = 0, 1, or 2 for a marker.

adm_estimates_test <- sapply(X = 1:ncol(alleles), FUN = algorithm_1snp_wrapper,
                        alleles = alleles, ancestries = adm_ancestries_test)

adm_estimates_test

