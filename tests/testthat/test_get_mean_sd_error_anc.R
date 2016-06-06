test_that("TESTING get_mean_sd_error_anc()", {

    print("Constructing an example where I know what to expect
          from this function.")

    estimates <- matrix(c(1:5,
                          0.15, 0.36, 0.57, 0.78, 0.99,
                          0.15, 0.36, 0.57, 0.78, 0.99,
                          0.15, 0.36, 0.57, 0.78, 0.99),
                        nrow = 5, ncol = 4)

    colnames(estimates) <- c("rsid", "Freq1_Anc0", "Freq1_Anc1", "Freq1_Anc2")
    print("Say estimates are: ")
    print(estimates)

    truth <- matrix(c(1:5,
                      0.1, 0.3, 0.5, 0.7, 0.9,
                      0.1, 0.3, 0.5, 0.7, 0.9,
                      0.1, 0.3, 0.5, 0.7, 0.9),
                    nrow = 5, ncol = 4)

    colnames(truth) <- c("rsid", "Freq1_Anc0", "Freq1_Anc1", "Freq1_Anc2")
    print("Say truth is: ")
    print(truth)

    errors <- estimates - truth

    # Test 1

    out1 <- get_mean_sd_error_anc(ancestry = 0,
                          estimates = estimates,
                          truth = truth)

    expect_equal(out1[1,], errors[,2])

    # Test 2

    out2 <- get_mean_sd_error_anc(ancestry = 1,
                                    estimates = estimates,
                                    truth = truth)

    expect_equal(out2[1,], errors[,2])

    # Test 3

    out3 <- get_mean_sd_error_anc(ancestry = 2,
                          estimates = estimates,
                          truth = truth)

    expect_equal(out3[1,], errors[,2])

})
