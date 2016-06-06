test_that("TESTING get_errors_summary_stats_1_scenario()", {

    print("The following should run: ")

    n_ind <- 3
    n_markers <- 2

    print(paste("n_ind = ", n_ind))
    print(paste("n_markers = ", n_markers))

    ancestries_matrix <- matrix(rep(c(0,0,1,1,2,2), times = 2),
                                nrow = 2 * n_ind, ncol = n_markers)

    print("ancestries_matrix = ")
    print(ancestries_matrix)

    # Rows: P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2).
    # Cols: SNPs.

    errors_matrix <- get_errors_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                                ancestries_matrix_true = ancestries_matrix,
                                                ancestries_matrix_estimated = ancestries_matrix)

    print("get_errors_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                                ancestries_matrix_true = ancestries_matrix,
                                                ancestries_matrix_estimated = ancestries_matrix)")
    print(errors_matrix)

    errors_summ_over_snps <- get_errors_summary_stats_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                        ancestries_matrix_true = ancestries_matrix,
                                        ancestries_matrix_estimated = ancestries_matrix)

    print("get_errors_summary_stats_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                        ancestries_matrix_true = ancestries_matrix,
                                        ancestries_matrix_estimated = ancestries_matrix)
          takes mean error over snp columns, and errors across snps should be identical,
          so mean errors should be the same as a column from the output of
          the get_errors_1_scenario() call")

    print(errors_summ_over_snps)

    expect_equal(errors_matrix[,1], errors_summ_over_snps[,1])
})

