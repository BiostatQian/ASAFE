# INPUTS:
#
# (1) Row of anc_spec_freqs:
# The particular [p0, p1, p2] setting we're considering.
#
# (2) anc_spec_freqs:
# Matrix with rows as [p0,p1,p2] settings,
# and columns (1) p0, (2) p1, (3) p2.
#
# (3) ancestries_matrix_true
# Rows: Individual chromosomes, 
# with chromosomes for the same individual consecutive.
# Cols: Markers.
# Entries: 0, 1, 2 for the true ancestry of origin.
#
# (4) ancestries_matrix_estimated
# Rows: Individual chromosomes, 
# with chromosomes for the same individual consecutive.
# Cols: Markers.
# Entries: 0, 1, 2 for the ancestry ASAFE will use to get estimates with.
# These ancestries may not be true ancestries.
#
# FUNCTION:
# Get errors for each ancestry-specific allele frequencies scenario.
#
# OUTPUTS:
#
# Matrix with
# Rows: [p0, p1, p2] settings.
# Cols:
# (1) p0, ancestry 0 specific allele 1 frequency
# (2) p1, ancestry 1 specific allele 1 frequency
# (3) p2, ancestry 2 specific allele 1 frequency
# (4) Mean of errors over SNPs in estimating P(Allele 1 | Anc 0)
# (5) SD of errors over SNPs in estimating P(Allele 1 | Anc 0)
# (6) Mean of errors over SNPs in estimating P(Allele 1 | Anc 1)
# (7) SD of errors over SNPs in estimating P(Allele 1 | Anc 1)
# (8) Mean of errors over SNPs in estimating P(Allele 1 | Anc 2)
# (9) SD of errors over SNPs in estimating P(Allele 1 | Anc 2)

get_scenario_errors <- function(row, anc_spec_freqs,  ancestries_matrix_true,
                                                      ancestries_matrix_estimated){

    # message("In get_scenario_errors(). 
    #          The first few elements of ancestries_matrix are: ")
    # message("ancestries_matrix[1,1:10]")
    # message(ancestries_matrix_estimated[1,])

    p0 <- anc_spec_freqs[row,1]
    p1 <- anc_spec_freqs[row,2]
    p2 <- anc_spec_freqs[row,3]

    # message(paste("row is ", row))
    # message(paste("p0 is ", p0))
    # message(paste("p1 is ", p1))
    # message(paste("p2 is ", p2))

    # 3 x 2 matrix.
    # Rows: (1) P(Allele 1 | Anc 0), (2) P(Allele 1 | Anc 1), (3) P(Allele 1 | Anc 2).
    # Cols: Mean errors over markers, SD of errors over markers.

    error_summary_stats <- suppressMessages(get_errors_summary_stats_1_scenario(p0, p1, p2,
                                                               ancestries_matrix_true,
                                                               ancestries_matrix_estimated))

    # 3 x 3 matrix: 3 rows of [p0, p1, p2].

    freqs_matrix <- matrix(c(p0, p1, p2),
                           nrow = nrow(error_summary_stats),
                           ncol = 3,
                           byrow = TRUE)

    colnames(freqs_matrix) <- c("p0", "p1", "p2")

    # cbind(3 x 3 matrix, 3 x 2 matrix), so a 3 x 5 matrix.

    out <- cbind(freqs_matrix, error_summary_stats)

    return(out)
}
