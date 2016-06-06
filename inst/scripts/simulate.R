#!/usr/bin/env Rscript
#
# Qian Zhang
# March 8, 2016
# 
# Simulate ancestries and genotypes.
#
# n_ind <- 250
# n_markers <- 50
#
# Call to run this script: Rscript 1 250 1000

# Prevent conflicts
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
jobID <- as.numeric(args[1])
n_ind <- as.numeric(args[2])
n_markers <- as.numeric(args[3])

print(date())
print(sessionInfo())

# Install R packages

# install.packages("../ASAFE_0.1.0.tar.gz", repos = NULL, type="source")

# This is what I ran on 3/11/2016
# library("ASAFE") 

# This is what I ran on 3/12/2016
library(ASAFE, lib.loc="/projects/geneva/gcc-fs2/R_packages/library")

# Ancestry-specific allele frequencies.
# Rows: Scenarios. 
# Cols: P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
 
anc_spec_freqs <- matrix(c(0, 0.1, 0.4, 0.48, 0.1, 0.4,
		         0.5, 0.5, 0.5, 0.5, 0.6, 0.5,
			 0.9, 0.8, 0.6, 0.52, 0.6, 0.5),
			 nrow = 6, ncol = 3)

#### Draw ancestries

set.seed(1)

# Ancestries matrix. 
# Rows: Individuals' alleles. Cols: Markers.

ancestries_matrix <- matrix(NA, nrow = 2 * n_ind, ncol = n_markers)
 
# Ancestries are 1, 2, or 3.

for(snp in 1:n_markers){
    
    ancestries_matrix[,snp] <- replicate(n = 2 * n_ind, sample(x = 0:2, size = 1))

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

ancestries_matrix <- ancestries_matrix[, col_all_3_ancs]

# In case I write over ancestries_matrix, I have the original
ancestries_matrix_original <- ancestries_matrix

#### Draw genotypes

# Draw an allele given ancestry and vector of ancestry-specific allele frequencies.

draw_allele_given_anc <- function(anc, freqs){
    
    if(anc == 0){
        allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[1], freqs[1]))
    }else if(anc == 1){
	allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[2], freqs[2]))
    }else if(anc == 2){
        allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[3], freqs[3]))
    }else{
 	stop('anc must be 0, 1, or 2')
    }

    return(allele)
}

##### Test

draw_allele_given_anc(anc = 1, freqs = anc_spec_freqs[1,])

# Alleles matrix. 
# Rows: Individuals' alleles. Cols: Markers.

alleles_matrix <- matrix(mapply(FUN = function(anc){
				  draw_allele_given_anc(anc, freqs = anc_spec_freqs[1,])
		               }, ancestries_matrix),
			 nrow = nrow(ancestries_matrix),
		         ncol = ncol(ancestries_matrix))

# Test that mapply works in the way that I think.
# 
# matrix(mapply(FUN = function(anc){
#                 anc + 1 # draw_allele_given_anc(anc, freqs = anc_spec_freqs[1,])
#              }, ancestries_matrix),
#       nrow = nrow(ancestries_matrix),
#       ncol = ncol(ancestries_matrix))

# To every locus, apply ASAFE.
estimates_1snp <- algorithm_1snp(alleles_1 = alleles_matrix[,1], 
		         	 ancestries_1 = ancestries_matrix[,1])

# Get estimated allele frequencies for all markers
# Rows: Ancestries 0, 1, or 2. Cols: Markers.

estimates <- sapply(X = 1:ncol(alleles_matrix),
		    FUN = function(marker_col){
		        algorithm_1snp(alleles_1 = alleles_matrix[,marker_col],
                                       ancestries_1 = ancestries_matrix[,marker_col])
		    })

#### END Test

# Get true ancestry-specific allele frequencies

get_true_freqs_1snp <- function(alleles_1, ancestries_1){
    
    indices_anc_0 <- which(ancestries_1 == 0)
    anc0_freq <- sum(alleles_1[indices_anc_0]) / length(indices_anc_0)
 
    indices_anc_1 <- which(ancestries_1 == 1)
    anc1_freq <- sum(alleles_1[indices_anc_1]) / length(indices_anc_1)

    indices_anc_2 <- which(ancestries_1 == 2)
    anc2_freq <- sum(alleles_1[indices_anc_2]) / length(indices_anc_2)   

    return(c(anc0_freq, anc1_freq, anc2_freq)) 
}

##### Test

true_freqs_1snp <- get_true_freqs_1snp(alleles_1 = alleles_matrix[,1], 
	          	   		  ancestries_1 = ancestries_matrix[,1])

# Get true allele frequencies for all markers.
# Rows: Ancestries 0, 1, or 2. Cols: Markers.

true_freqs <- sapply(X = 1:ncol(alleles_matrix),
                     FUN = function(marker_col){
                        get_true_freqs_1snp(alleles_1 = alleles_matrix[,marker_col],
                                            ancestries_1 = ancestries_matrix[,marker_col])
                     })  

# Get errors.
# Rows: Ancestries 0, 1, or 2. Cols: Markers.

errors <- estimates - true_freqs

# Test that matrix subtraction works element-wise in R:
matrix(2:5, nrow = 2, ncol = 2) - matrix(1, nrow = 2, ncol = 2)

# Mean (across markers) error for ancestries 0, 1, and 2:

mean_errors <- apply(X = errors, MARGIN = 1, FUN = mean)
sd_errors <- apply(X = errors, MARGIN = 1, FUN = sd)

##### End Test

# Given a vector of ancestries and a vector of alleles, both listed in this order:
# [ID 1 Allele 1, ID 1 Allele 2, ..., ID n Allele 1, ID n Allele 2],
# count the number of instances where the individual was both het ancestry and het genotype.

count_het_anc_het_geno <- function(alleles_1, ancestries_1){

   count_hets <- 0
   
   num_ind <- length(alleles_1) / 2

   for(ind in 1:num_ind){
       
       anc_pair <- ancestries_1[c(2*ind - 1, 2*ind)]
       genotype <- alleles_1[c(2*ind - 1, 2*ind)]

       if(length(unique(anc_pair)) == 2 && 
	  length(unique(genotype)) == 2){

           count_hets <- count_hets + 1

       }

   }

   return(count_hets)    
}

###### Test 

count_het_anc_het_geno(alleles_1 = alleles_matrix[,1], ancestries_1 = ancestries_matrix[,1])

###### END Test

# Make a function that will draw alleles given ancestries, and then
# run ASAFE on the alleles and ancestries data, and then get
# Rows: P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA), Number of individuals with het anc and geno.
# Cols: Markers
# for a particular set of ancestral allele frequencies p0, p1, and p2.

get_errors_1_scenario <- function(p0, p1, p2, ancestries_matrix_true,
                                                ancestries_matrix_estimated){

    # Alleles matrix.
    # Rows: Individuals' alleles. Cols: Markers.

    alleles_matrix <- matrix(mapply(FUN = function(anc){
                                        draw_allele_given_anc(anc, freqs = c(p0, p1, p2))
                                    }, ancestries_matrix_true),
                                    nrow = nrow(ancestries_matrix_true),
                                    ncol = ncol(ancestries_matrix_true))

    # Number of individuals with both het ancestry and het genotype, for all markers.
    # Vector of length the number of markers.
    het_counts <- sapply(X = 1:ncol(alleles_matrix),
                         FUN = function(marker_col){
                             count_het_anc_het_geno(alleles_1 = alleles_matrix[,marker_col],
                                                    ancestries_1 = ancestries_matrix_true[,marker_col])
                         })

    # 3 x Number of marker matrix of ancestry-specific allele frequencies per marker
    estimates <- sapply(X = 1:ncol(alleles_matrix),
                        FUN = function(marker_col){
                            algorithm_1snp(alleles_1 = alleles_matrix[,marker_col],
                                           ancestries_1 = ancestries_matrix_estimated[,marker_col])
                        })

    # 3 x Number of marker matrix of true ancestry-specific allele frequencies.
    true_freqs <- sapply(X = 1:ncol(alleles_matrix),
                         FUN = function(marker_col){
                             get_true_freqs_1snp(alleles_1 = alleles_matrix[,marker_col],
                                                 ancestries_1 = ancestries_matrix_true[,marker_col])
                         })

    # 3 x (Number of markers) matrix
    errors <- (estimates - true_freqs)

    # 3 x (Number of markers) matrix
    relative_errors <- (estimates - true_freqs) / (true_freqs) * 100

    # Rows: P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA), Number of individuals with het anc and geno
    # Cols: Markers

    out <- rbind(errors, relative_errors, het_counts)

    return(out)
}

##### Test

errors_counts_matrix <- get_errors_1_scenario(p0 = 0.1, p1 = 0.5, p2 = 0.9, 
						ancestries_matrix_true = ancestries_matrix_original,
						ancestries_matrix_estimated = ancestries_matrix_original)

##### End TEst

# Inputs:
# (1)-(3) Ancestral P(Allele 1 | Anc) frequencies, for Anc = 0, 1, or 2.
# (4) ancestries_matrix_true: Matrix of true ancestries
# (rows are alleles with alleles for the same individual consecutive, cols are markers)
# (5) ancestries_matrix_eximated: Matrix of estimated ancestries
# that ASAFE will use
#
# Outputs:
# (1) A 7 x 2 matrix with:
# - Rows: Errors in P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA), and then
# Relative Errors in those ancestry-specific allele frequencies, 
# Number of individuals with het ancestry and het genotype
# - Cols: Mean and SD across SNPs

get_errors_summary_stats_1_scenario <- function(p0, p1, p2, ancestries_matrix_true,
                                                            ancestries_matrix_estimated){

    # message("In get_errors_summary_stats_1_scenario().")
    # message("ancestries_matrix[1,] is: ")
    # message(ancestries_matrix[1,])

    errors_het_counts <- get_errors_1_scenario(p0, p1, p2, ancestries_matrix_true,
                                                            ancestries_matrix_estimated)

    # Mean (across markers) error for ancestries 0, 1, and 2:

    mean_errors_het_counts <- apply(X = errors_het_counts, MARGIN = 1, FUN = mean)
    sd_errors_het_counts <- apply(X = errors_het_counts, MARGIN = 1, FUN = sd)

    output <- cbind(mean_errors_het_counts, sd_errors_het_counts)

    return(output)
}

##### Test

get_errors_summary_stats_1_scenario(p0 = 0.5, p1 = 0.5, p2 = 0.5, 
					ancestries_matrix_true = ancestries_matrix_original,
					ancestries_matrix_estimated = ancestries_matrix_original)

##### END Test

# Get errors for each ancestry-specific allele frequencies scenario.
get_scenario_errors <- function(row, anc_spec_freqs,  ancestries_matrix_true,
                                                      ancestries_matrix_estimated){

    message("In get_scenario_errors(). The first few elements of ancestries_matrix are: ")
    message("ancestries_matrix[1,1:10]")
    message(ancestries_matrix_estimated[1,])

    p0 <- anc_spec_freqs[row,1]
    p1 <- anc_spec_freqs[row,2]
    p2 <- anc_spec_freqs[row,3]

    message(paste("row is ", row))
    message(paste("p0 is ", p0))
    message(paste("p1 is ", p1))
    message(paste("p2 is ", p2))

    error_summary_stats <- suppressMessages(get_errors_summary_stats_1_scenario(p0, p1, p2,
                                                               ancestries_matrix_true,
                                                               ancestries_matrix_estimated))

    freqs_matrix <- matrix(c(p0, p1, p2),
                           nrow = nrow(error_summary_stats),
                           ncol = 3,
                           byrow = TRUE)

    colnames(freqs_matrix) <- c("p0", "p1", "p2")

    out <- cbind(freqs_matrix, error_summary_stats)

    return(out)
}

##### Test

scenario_1 <- get_scenario_errors(row = 1, anc_spec_freqs = anc_spec_freqs,
				  ancestries_matrix_true = ancestries_matrix_original, 
				  ancestries_matrix_estimated = ancestries_matrix_original)
scenario_1

scenario_2 <- get_scenario_errors(row = 2, anc_spec_freqs = anc_spec_freqs,
                                  ancestries_matrix_true = ancestries_matrix_original, 
                                  ancestries_matrix_estimated = ancestries_matrix_original)
scenario_2

get_scenario_errors(row = 3, anc_spec_freqs = anc_spec_freqs,
                                  ancestries_matrix_true = ancestries_matrix_original, 
                                  ancestries_matrix_estimated = ancestries_matrix_original)

get_scenario_errors(row = 4, anc_spec_freqs = anc_spec_freqs,
                                  ancestries_matrix_true = ancestries_matrix_original, 
                                  ancestries_matrix_estimated = ancestries_matrix_original)

get_scenario_errors(row = 5, anc_spec_freqs = anc_spec_freqs,
                                  ancestries_matrix_true = ancestries_matrix_original, 
                                  ancestries_matrix_estimated = ancestries_matrix_original)

get_scenario_errors(row = 6, anc_spec_freqs = anc_spec_freqs,
                                  ancestries_matrix_true = ancestries_matrix_original, 
                                  ancestries_matrix_estimated = ancestries_matrix_original)

#### END Test

###### Add in ancestry-calling error.
# Change each ancestry to a different ancestry with equal probability,
# changing to the other thing with probability 7%.

sample_ancestry <- function(anc){
    
    if(anc == 0){ 
        out <- sample(x = c(1,2), size = 1, prob = c(0.5, 0.5))
    }else if(anc == 1){
	out <- sample(x = c(0,2), size = 1, prob = c(0.5, 0.5)) 
    }else if(anc == 2){
	out <- sample(x = c(0,1), size = 1, prob = c(0.5, 0.5))
    }else{
	stop("anc should be 0, 1, or 2")
    }

    return(out)
}

#### Test

sample_ancestry(anc = 0)

#### End Test

# Change ancestry. Decide to change ancestry with probability error_rate.
# If decide to change ancestry, then sample the other two ancestries with 50/50 chance.

change_ancestry <- function(anc, error_rate){
    
    ind_change <- sample(x = 0:1, size = 1, prob = c(1 - error_rate, error_rate))  
   
    if(ind_change == 1){
	message("ind_change == 1")
        out <- sample_ancestry(anc)
    }else{
	out <- anc
    }
 
    return(out)
}

#### Test

change_ancestry(anc = 0, error_rate = 0.07)

#### End Test

#### Get accuracy with the ancestries matrix that has errors in it.

# Inputs:
# (1) Error_rate: P(We change a true ancestry value)
# (2) anc_spec_freqs: Scenarios x 3 ancestry-specific allele frequencies matrix
# (3) ancestries_matrix_true: Alleles i.e. 2 * number of individuals x markers matrix of true ancestries (0,1,or 2)
# 
# Outputs:
# (1) For each combo of true ancestral allele frequencies,
# give 7 summary stats (rows) with cols (p0, p1, p2, mean, sd).

get_results_error <- function(error_rate, anc_spec_freqs, ancestries_matrix_true){

    # Introduce ancestry errors.

    ancestries_error_matrix <- matrix(mapply(FUN = function(anc){
                                                       suppressMessages(change_ancestry(anc, error_rate))
                                                   }, ancestries_matrix_true),
                                                   nrow = nrow(ancestries_matrix_true),
                                                   ncol = ncol(ancestries_matrix_true))

    out <- matrix(NA, nrow = 7 * nrow(anc_spec_freqs), ncol = 5)

    # For each scenario, get the table of error summary stats and double het counts
    # 4 rows per scenarios: P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA), het count
    # Columns: P(Allele 1 in Ancestral Pop 1), P(Allele 1 in Ancestral Pop 2),
    #          P(Allele 1 in Ancestral Pop 3), Mean, SD.

    for(row in 1:nrow(anc_spec_freqs)){

        summary_stats <- get_scenario_errors(row, anc_spec_freqs,
                                             ancestries_matrix_true = ancestries_matrix_true,
                                             ancestries_matrix_estimated = ancestries_error_matrix)

        out[(row - 1) * 7 + 1:7, ] <- summary_stats

    }

    colnames(out) <- c("p0", "p1", "p2", "Mean", "SD")

    row.names(out) <- replicate(n = nrow(anc_spec_freqs),
				          expr = c("Abs_Error_Afr_Freq",          
                                                   "Abs_Error_Eur_Freq",          
                                                   "Abs_Error_NA_Freq",         
                                                   "Rel_Error_Afr_Freq",
                                                   "Rel_Error_Eur_Freq",
                                                   "Rel_Error_NA_Freq",
                                                   "Number_Ind_Double_Het"))

    return(out)
}

# Error rate = 0

error_rate <- 0

out_0.0 <- get_results_error(error_rate = error_rate, 
                             anc_spec_freqs = anc_spec_freqs,
                             ancestries_matrix = ancestries_matrix_original)

write.table(x = out_0.0, file = paste("error_summary_stats_error_", error_rate,
			 "_jobID_", jobID,	 
			 ".txt", sep = ""),
                         quote = FALSE, 
                         col.names = TRUE, row.names = TRUE, append = FALSE)

# Error rate = 0.07

error_rate <- 0.07

out_0.07 <- get_results_error(error_rate = error_rate, 
                         anc_spec_freqs = anc_spec_freqs,
                         ancestries_matrix = ancestries_matrix_original) 

write.table(x = out_0.07, file = paste("error_summary_stats_error_", error_rate,
                         "_jobID_", jobID,
                         ".txt", sep = ""),
                         quote = FALSE,
                         col.names = TRUE, row.names = TRUE, append = FALSE)

# Error rate = 1

error_rate <- 1.0

out_1.0 <- get_results_error(error_rate = error_rate, 
                             anc_spec_freqs = anc_spec_freqs,
                             ancestries_matrix = ancestries_matrix_original)

write.table(x = out_1.0, file = paste("error_summary_stats_error_", error_rate,
                         "_jobID_", jobID,
                         ".txt", sep = ""),
                         quote = FALSE,
                         col.names = TRUE, row.names = TRUE, append = FALSE)

print(date())
