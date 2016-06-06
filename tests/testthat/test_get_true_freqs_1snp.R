test_that("TESTING get_true_freqs_1snp()", {

    print("get_true_freqs_1snp(alleles_1 = c(0, 0, 1, 0, 1, 1),
           ancestries_1 = c(0, 1, 0, 2, 1, 2))
           should return [0.5, 0.5, 0.5]")

    expect_equal(get_true_freqs_1snp(alleles_1 = c(0, 0, 1, 0, 1, 1),
                        ancestries_1 = c(0, 1, 0, 2, 1, 2)),
                 c(0.5, 0.5, 0.5))

})
