# Inputs:
#
# (1)-(3) Ancestral P(Allele 1 | Anc) frequencies, for Anc = 0, 1, or 2.
#
# (4) ancestries_matrix_true: 
# Matrix of true ancestries
# (rows are alleles with alleles for the same individual consecutive, 
# cols are markers)
#
# (5) ancestries_matrix_eximated: 
# Matrix of estimated ancestries that ASAFE will use
# (rows are alleles with alleles for the same individual consecutive, 
# cols are markers)
#
# Outputs:
# (1) A 3 x 2 matrix with:
# - Rows: Errors in P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA)
# - Cols: Mean and SD of errors across SNPs

get_errors_summary_stats_1_scenario <- function(p0, p1, p2, 
                                                ancestries_matrix_true,
                                                ancestries_matrix_estimated){

    # Matrix with:
    # Rows: P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
    # Cols: SNPS.
    errors <- get_errors_1_scenario(p0, p1, p2, ancestries_matrix_true,
                                                            ancestries_matrix_estimated)

    # Mean (across markers) error for ancestries 0, 1, and 2:
    
    mean_errors <- apply(X = errors, MARGIN = 1, FUN = mean)
    sd_errors <- apply(X = errors, MARGIN = 1, FUN = sd)

    # 3 x 2 matrix: 
    # Rows: Things we're trying to estimate,
    # P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
    # Cols: Mean error over SNPs, SD of errors over SNPs.

    output <- cbind(mean_errors, sd_errors)

    return(output)
}
