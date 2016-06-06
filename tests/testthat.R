# The line you need to have R CMD check work,
# and avoid some weird error I got from R CMD check.

Sys.setenv("R_TESTS" = "")

# The following commands model off of 
# https://github.com/hadley/testthat

library(testthat)
library(ASAFE)

test_check("ASAFE")

# Note that in the directory tests/testthat,
# there's 1 test file per function. I like that format,
# because I don't have to scroll through a long file.
# Hadley does the same thing here: 
# https://github.com/hadley/testthat/tree/master/tests/testthat.
