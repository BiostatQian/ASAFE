test_that("TESTING get_scenario_errors()", {

    anc_spec_freqs <- matrix(c(1, 1, 1,
                               0, 0, 0,
                               0.5, 0.5, 0.5), nrow = 3, ncol = 3, byrow = TRUE)

    print("anc_spec_freqs = ")
    print(anc_spec_freqs)

    ancestries_matrix <- matrix(c(0, 0, 1, 1, 2, 2),
                                nrow = 6, ncol = 2)

    print("ancestries_matrix = ")
    print(ancestries_matrix)

    out1 <- get_scenario_errors(row = 1, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)

    print("get_scenario_errors(row = 1, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)")
    print(out1)

    out2 <- get_scenario_errors(row = 2, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)

    print("get_scenario_errors(row = 2, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)")

    print(out2)

    out3 <- get_scenario_errors(row = 3, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)

    print("get_scenario_errors(row = 3, anc_spec_freqs = anc_spec_freqs,
                        ancestries_matrix_true = ancestries_matrix,
                        ancestries_matrix_estimated = ancestries_matrix)")

    print(out3)
})

