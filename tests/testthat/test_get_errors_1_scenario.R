test_that("TESTING get_errors_1_scenario()", {

    print("The following should run: ")

    n_ind <- 3
    n_markers <- 2

    print(paste("n_ind = ", n_ind))
    print(paste("n_markers = ", n_markers))

    ancestries_matrix <- matrix(rep(c(0,0,1,1,2,2), times = 2),
                                nrow = 2 * n_ind, ncol = n_markers)

    print("ancestries_matrix = ")
    print(ancestries_matrix)

    errors_counts_matrix <- get_errors_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                                ancestries_matrix_true = ancestries_matrix,
                                                ancestries_matrix_estimated = ancestries_matrix)

    print("get_errors_1_scenario(p0 = 1, p1 = 1, p2 = 1,
                                                ancestries_matrix_true = ancestries_matrix,
                                                ancestries_matrix_estimated = ancestries_matrix)")
    print(errors_counts_matrix)

})

