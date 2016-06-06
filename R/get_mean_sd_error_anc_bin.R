# get_mean_sd_error_anc_bin: Get mean error and sd error
# for a particular ancestry's allele 1 frequency estimates,
# for SNPs with true allele 1 frequency in a particular bin (min_bin, max_bin)
#
# Inputs:
# (1) Ancestry (0, 1, or 2). This is a scalar.
# (2) Min allele frequency min_bin of bin. This is a scalar.
# (3) Max allele frequency max_bin of bin. This is a scalar.
# (4) Estimated allele freqencies. This is a matrix.
# - Rows: Markers.
# - Cols: Marker, and then
# estimated P(A1 | Anc 0), P(A1 | Anc 1), and P(A1 | Anc 2).
# (5) True allele frequencies. This is a matrix.
# - Rows: Markers.
# - Cols: Marker, and then
# true P(A1 | Anc 0), P(A1 | Anc 1), and P(A1 | Anc 2).
#
# Outputs:
# (1) Mean
# (2) SD
# of the errors (where error = estimate - truth) for markers with
# true ancestry specific allele frequency in bin (min_bin, max_bin].

get_mean_sd_error_anc_bin <- function(ancestry, min_bin, max_bin,
                                      estimates, truth){

    errors <- estimates[,ancestry+2] - truth[,ancestry+2]

    indices <- intersect(which(truth[,ancestry+2] > min_bin),
                        which(truth[,ancestry+2] <= max_bin))

    errors_bin <- errors[indices]

    out <- c(mean(errors_bin), sd(errors_bin))

    return(out)

}
