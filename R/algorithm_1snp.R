# Function to apply EM algorithm to 1 SNP
#
# INPUTS
#
# 1) alleles_1: 
# Vector of alleles (0 or 1) for each individual's 2 chromosomes,
# with chromosomes for the same individual consecutive.
#
# Example: Alleles in order ADM1, ADM1, ADM2, ADM2, ..., ADM250, ADM250.
#
# 2) ancestries_1: 
# Vector of ancestries (0, 1, or 2) for each individual's 2 chromosomes,
# with chromosomes for the same individual consecutive.
#
# Example: Ancestries in order ADM1, ADM1, ADM2, ADM2, ..., ADM250, ADM250.
#
# OUTPUTS
#
# 1) Ancestry-specific allele frequency estimates from the EM Algorithm:
# P(1|A), P(1|E), P(1|N).
#
# DETAILS:
# Observed data vector n is called mprime in the ASAFE paper supplement.

algorithm_1snp <- function(alleles_1, ancestries_1){

    # Create observed data category counts n based on
    # ancestry and genotype data.

    # Initialize counts.
    n <- rep(0, times = 18)

    # Update counts, iterating over individuals.
    for(i in 1:(length(alleles_1)/2)){
        ancestries_ind <- c(ancestries_1[2*i-1], ancestries_1[2*i])
        alleles_ind <- c(alleles_1[2*i - 1], alleles_1[2*i])

        n <- increment_n(n, anc_ind = ancestries_ind, geno_ind = alleles_ind)
    }

    ######## Apply EM to observed counts

    p1_adm <- em(n = n, epsilon = 10^-8, iteration_cap = 1000)

    return(p1_adm)
}
