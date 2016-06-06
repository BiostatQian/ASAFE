test_that("TESTING get_results_error()", {

    n_ind <- 3
    n_markers <- 2

    print(paste("n_ind = ", n_ind))
    print(paste("n_markers = ", n_markers))

    anc_spec_freqs <- matrix(c(1, 1, 1,
                               0, 0, 0,
                               0.5, 0.5, 0.5), nrow = 3, ncol = 3, byrow = TRUE)

    print("anc_spec_freqs = ")
    print(anc_spec_freqs)

    ancestries_matrix <- matrix(rep(c(0,0,1,1,2,2), times = 2),
                                nrow = 2 * n_ind, ncol = n_markers)

    print("ancestries_matrix = ")
    print(ancestries_matrix)

    print("get_results_error(error_rate = 0,
           anc_spec_freqs, ancestries_matrix_true) = ")

    print(get_results_error(error_rate = 0,
                      anc_spec_freqs,
                      ancestries_matrix_true = ancestries_matrix))

    print("get_results_error(error_rate = 1,
           anc_spec_freqs, ancestries_matrix_true) = ")

    print(get_results_error(error_rate = 1,
                      anc_spec_freqs,
                      ancestries_matrix_true = ancestries_matrix))

})

