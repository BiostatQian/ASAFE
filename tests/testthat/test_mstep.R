test_that("TESTING mstep(): M-Step returns a numeric vector 6 elements long", {

    expect_equal(length(mstep(m = rep(1, times = 21))), 6)

    expect_match(class(mstep(m = rep(1, times = 21))), "numeric")

})

