#### Get accuracy with the ancestries matrix that has errors in it.
#
# Inputs:
#
# (1) Error_rate: P(We change a true ancestry value)
#
# (2) anc_spec_freqs: 
# [p0, p1, p2] scenarios x 3 ancestry-specific allele frequencies matrix
#
# (3) ancestries_matrix_true: 
# 1, ..., 2*(number of individuals) alleles 
# x 1, ..., (number of markers) matrix of true ancestries (0,1,or 2)
# 
# Outputs:
#
# (1) For each combo of true ancestral allele frequencies [p0, p1, p2],
# give summary stats (mean and SD) of errors over markers 
# in estimating each of 
# P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
#
# Output is a matrix with:
# Rows: P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2). 
# Cols: (p0, p1, p2, mean of errors over markers, sd of errors over markers).

get_results_error <- function(error_rate, 
                              anc_spec_freqs, 
                              ancestries_matrix_true){

    # Introduce ancestry errors.

    ancestries_error_matrix <- matrix(mapply(FUN = function(anc){
                                                       suppressMessages(change_ancestry(anc, error_rate))
                                                   }, ancestries_matrix_true),
                                      nrow = nrow(ancestries_matrix_true),
                                      ncol = ncol(ancestries_matrix_true))

    out <- matrix(NA, nrow = 3 * nrow(anc_spec_freqs), ncol = 5)

    # For each [p0, p1, p2] scenario, get the table of error summary stats.
    # 3 rows per scenario: 
    # P(Allele 1 | Afr), P(Allele 1 | Eur), P(Allele 1 | NA).
    # Columns: P(Allele 1 in Ancestral Pop 0), P(Allele 1 in Ancestral Pop 1),
    #          P(Allele 1 in Ancestral Pop 2), Mean, SD.

    for(row in 1:nrow(anc_spec_freqs)){

        # 3 x 5 matrix 
        summary_stats <- get_scenario_errors(row = row, 
                                             anc_spec_freqs = anc_spec_freqs,
                                             ancestries_matrix_true = ancestries_matrix_true,
                                             ancestries_matrix_estimated = ancestries_error_matrix)

        out[(row - 1) * 3 + 1:3, ] <- summary_stats

    }

    colnames(out) <- c("p0", "p1", "p2", "Mean", "SD")

    row.names(out) <- replicate(n = nrow(anc_spec_freqs),
                                          expr = c("Abs_Error_Afr_Freq",
                                                   "Abs_Error_Eur_Freq",
                                                   "Abs_Error_NA_Freq"))

    return(out)
}

