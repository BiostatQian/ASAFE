test_that("TESTING estep(): EStep returns a numeric vector with 21 entries", {

    expect_match(class(estep(mprime = rep(1, times = 18), p = rep(1/6, times = 6))),
                 "numeric")

    expect_equal(length(estep(mprime = rep(1, times = 18), p = rep(1/6, times = 6))),
                 21)

})


