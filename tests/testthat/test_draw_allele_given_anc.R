test_that("TESTING draw_allele_given_anc()", {

    # Should return 1
    print("draw_allele_given_anc(anc = 0, freqs = c(1, 0, 0)) should return 1")
    expect_equal(draw_allele_given_anc(anc = 0, freqs = c(1, 0, 0)), 1)

    # Should return 0
    print("draw_allele_given_anc(anc = 0, freqs = c(0, 0, 0)) should return 0")
    expect_equal(draw_allele_given_anc(anc = 0, freqs = c(0, 0, 0)), 0)

    # Should return 1
    print("draw_allele_given_anc(anc = 1, freqs = c(0, 1, 0)) should return 1")
    expect_equal(draw_allele_given_anc(anc = 1, freqs = c(0, 1, 0)), 1)

    # Should return 0
    print("draw_allele_given_anc(anc = 1, freqs = c(0, 0, 0)) should return 0")
    expect_equal(draw_allele_given_anc(anc = 1, freqs = c(0, 0, 0)), 0)

    # Should return 1
    print("draw_allele_given_anc(anc = 2, freqs = c(0, 0, 1)) should return 1")
    expect_equal(draw_allele_given_anc(anc = 2, freqs = c(0, 0, 1)), 1)

    # Should return 0
    print("draw_allele_given_anc(anc = 2, freqs = c(0, 0, 0)) should return 0")
    expect_equal(draw_allele_given_anc(anc = 2, freqs = c(0, 0, 0)), 0)

    # Should return "anc must be 0, 1, or 2"
    print("draw_allele_given_anc(anc = 3, freqs = c(0, 0, 0)) should return error")
    expect_error(draw_allele_given_anc(anc = 3, freqs = c(0, 0, 0)))

})

