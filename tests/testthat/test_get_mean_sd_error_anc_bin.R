test_that("TESTING get_mean_sd_error_anc_bin()", {

    print("Constructing an example where I know what to expect
          from this function.")

    estimates <- matrix(c(1:6,
                          0.1, 0.1, 0.1, 0.9, 0.9, 0.9,
                          0.1, 0.1, 0.1, 0.9, 0.9, 0.9,
                          0.1, 0.1, 0.1, 0.9, 0.9, 0.9),
                          nrow = 6, ncol = 4)

    colnames(estimates) <- c("rsid", "Freq1_Anc0", "Freq1_Anc1", "Freq1_Anc2")
    print("Say estimates are: ")
    print(estimates)

    truth <- matrix(c(1:6,
                    0.2, 0.2, 0.5, 0.5, 0.8, 0.8,
                    0.2, 0.2, 0.5, 0.5, 0.8, 0.8,
                    0.2, 0.2, 0.5, 0.5, 0.8, 0.8),
                    nrow = 6, ncol = 4)

    colnames(truth) <- c("rsid", "Freq1_Anc0", "Freq1_Anc1", "Freq1_Anc2")
    print("Say truth is: ")
    print(truth)

    errors <- estimates - truth

    # Test 1

    out1 <- get_mean_sd_error_anc_bin(ancestry = 0,
                                    min_bin = 0.1, max_bin = 0.3,
                                    estimates = estimates,
                                    truth = truth)

    expect_equal(out1[1], -0.1)

    # Test 2

    out2 <- get_mean_sd_error_anc_bin(ancestry = 1,
                                      min_bin = 0.3, max_bin = 0.6,
                                      estimates = estimates,
                                      truth = truth)
    expect_equal(out2[1], 0)

    # Test 3

    out3 <- get_mean_sd_error_anc_bin(ancestry = 1,
                                      min_bin = 0.7, max_bin = 0.9,
                                      estimates = estimates,
                                      truth = truth)
    expect_equal(out3[1], 0.1)

})
