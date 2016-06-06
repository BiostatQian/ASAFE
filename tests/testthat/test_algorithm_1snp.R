test_that("TESTING algorithm_1snp() ", {

    # Define alleles_1 and ancestries_1 such that
    # the counts of the numbers of observations
    # in each of 18 categories is a vector of 1's.
    # This should yield the same result as
    # em(n = replicate(n = 18, expr = 1),
    #          epsilon = 10^-8, iteration_cap = 1000).

    alleles_1_00 <- as.vector(replicate(n = 6, expr = c(0,0)))
    alleles_1_01 <- as.vector(replicate(n = 6, expr = c(0,1)))
    alleles_1_11 <- as.vector(replicate(n = 6, expr = c(1,1)))

    alleles_1 <- c(alleles_1_00, alleles_1_01, alleles_1_11)

    ancestries_1_part <- c(0,0, 0,1, 0,2, 1,1, 1,2, 2,2)

    ancestries_1 <- c(ancestries_1_part, ancestries_1_part, ancestries_1_part)

    out <- algorithm_1snp(alleles_1 = alleles_1,
                          ancestries_1 = ancestries_1)

    print("alleles_1 and ancestries_1 have been defined so that
          n is an 18-long vector of 1's,
          algorithm_1snp(alleles_1, ancestries_1) =
          em(n = replicate(n = 18, expr = 1), epsilon = 10^-8, iteration_cap = 1000).
          algorithm_1snp(alleles_1, ancestries_1) = ")

    print(out)

    expect_equal(out, c(0.5, 0.5, 0.5))
})

