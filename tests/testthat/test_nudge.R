test_that("TESTING nudge: nudge() only moves probabilities away from 0 or 1", {

    expect_equal(nudge(p = 0), 1e-08)
    expect_equal(nudge(p = 1), 1 - 10^-8)
    expect_equal(nudge(p = 0.5), 0.5)

})
