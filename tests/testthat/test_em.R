test_that("TESTING em()", {

    print("em(n = replicate(n = 18, expr = 1),
       epsilon = 10^-8, iteration_cap = 1000) = ")

    out <- em(n = replicate(n = 18, expr = 1),
              epsilon = 10^-8, iteration_cap = 1000)

    print(out)

    expect_equal(length(out), 3)

})
