test_that("TESTING change_ancestry()", {

    print("change_ancestry(anc = 0, error_rate = 1) should be 1 or 2")
    print(change_ancestry(anc = 0, error_rate = 1))

    expect_equal(change_ancestry(anc = 0, error_rate = 0), 0)

    print("change_ancestry(anc = 1, error_rate = 1) should be 0 or 2")
    print(change_ancestry(anc = 1, error_rate = 1))

    expect_equal(change_ancestry(anc = 1, error_rate = 0), 1)

    print("change_ancestry(anc = 2, error_rate = 1) should be 0 or 1")
    print(change_ancestry(anc = 2, error_rate = 1))

    expect_equal(change_ancestry(anc = 2, error_rate = 0), 2)
})
