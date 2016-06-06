# get_mean_sd_error_anc:
#
# INPUTS:
# (1) ancestry (0,1, or 2)
# (2) estimates: Matrix with rows as markers and
# columns as marker ID, estimated ancestry-specific allele 1 frequencies
# for ancestries 0, 1, and 2.
# (3) truth: Matrix with rows as markers and
# columns as marker ID, true ancestry-specific allele 1 frequencies
# for ancestries 0, 1, and 2.
#
# OUTPUTS:
# (1) 2 rows like in Table 1 of ASAFE paper for a certain ancestry.
# Rows are mean error, sd error.
# Cols are allele frequency bins.

get_mean_sd_error_anc <- function(ancestry, estimates, truth){

    out <- matrix(nrow = 2, ncol = 5)

    out[,1] <- get_mean_sd_error_anc_bin(ancestry = ancestry,
                                         min_bin = 0, max_bin = 0.2,
                                         estimates = estimates, truth = truth)

    out[,2] <- get_mean_sd_error_anc_bin(ancestry = ancestry,
                                         min_bin = 0.2, max_bin = 0.4,
                                         estimates = estimates, truth = truth)

    out[,3] <- get_mean_sd_error_anc_bin(ancestry = ancestry,
                                         min_bin = 0.4, max_bin = 0.6,
                                         estimates = estimates, truth = truth)

    out[,4] <- get_mean_sd_error_anc_bin(ancestry = ancestry,
                                         min_bin = 0.6, max_bin = 0.8,
                                         estimates = estimates, truth = truth)

    out[,5] <- get_mean_sd_error_anc_bin(ancestry = ancestry,
                                         min_bin = 0.8, max_bin = 1,
                                         estimates = estimates, truth = truth)

    return(out)
}
