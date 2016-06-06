test_that("TESTING increment_n: See which element of vector n gets changed by giving increment_n
                 a certain genotype or ancestry pair", {

    # Test genotype 0/0
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,0), geno_ind = c(0, 0)) == 1), 1)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,1), geno_ind = c(0, 0)) == 1), 4)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,2), geno_ind = c(0, 0)) == 1), 6)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,0), geno_ind = c(0, 0)) == 1), 2)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,1), geno_ind = c(0, 0)) == 1), 2)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,0), geno_ind = c(0, 0)) == 1), 3)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,2), geno_ind = c(0, 0)) == 1), 3)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,1), geno_ind = c(0, 0)) == 1), 5)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,2), geno_ind = c(0, 0)) == 1), 5)

    # Test genotype 1/0
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,0), geno_ind = c(1, 0)) == 1), 7)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,1), geno_ind = c(1, 0)) == 1), 10)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,2), geno_ind = c(1, 0)) == 1), 12)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,0), geno_ind = c(1, 0)) == 1), 8)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,1), geno_ind = c(1, 0)) == 1), 8)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,0), geno_ind = c(1, 0)) == 1), 9)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,2), geno_ind = c(1, 0)) == 1), 9)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,1), geno_ind = c(1, 0)) == 1), 11)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,2), geno_ind = c(1, 0)) == 1), 11)

    # Test genotype 0/1
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,0), geno_ind = c(0, 1)) == 1), 7)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,1), geno_ind = c(0, 1)) == 1), 10)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,2), geno_ind = c(0, 1)) == 1), 12)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,0), geno_ind = c(0, 1)) == 1), 8)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,1), geno_ind = c(0, 1)) == 1), 8)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,0), geno_ind = c(0, 1)) == 1), 9)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,2), geno_ind = c(0, 1)) == 1), 9)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,1), geno_ind = c(0, 1)) == 1), 11)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,2), geno_ind = c(0, 1)) == 1), 11)

    # Test genotype 1/1
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,0), geno_ind = c(1, 1)) == 1), 13)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,1), geno_ind = c(1, 1)) == 1), 16)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,2), geno_ind = c(1, 1)) == 1), 18)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,0), geno_ind = c(1, 1)) == 1), 14)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,1), geno_ind = c(1, 1)) == 1), 14)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,0), geno_ind = c(1, 1)) == 1), 15)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(0,2), geno_ind = c(1, 1)) == 1), 15)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(2,1), geno_ind = c(1, 1)) == 1), 17)
    expect_equal(which(increment_n(n = rep(0, times = 18), anc_ind = c(1,2), geno_ind = c(1, 1)) == 1), 17)

})
