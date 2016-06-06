test_that("TESTING algorithm_1snp_wrapper()", {

    alleles_1_00 <- as.vector(replicate(n = 6, expr = c(0,0)))
    alleles_1_01 <- as.vector(replicate(n = 6, expr = c(0,1)))
    alleles_1_11 <- as.vector(replicate(n = 6, expr = c(1,1)))
    alleles_1 <- c(alleles_1_00, alleles_1_01, alleles_1_11)

    ancestries_1_part <- c(0,0, 0,1, 0,2, 1,1, 1,2, 2,2)
    ancestries_1 <- c(ancestries_1_part, ancestries_1_part, ancestries_1_part)

    alleles <- cbind(alleles_1 , alleles_1)
    colnames(alleles) <- c("rs1", "rs2")

    ancestries <- rbind(ancestries_1, ancestries_1)

    print("alleles are: ")
    print(alleles)

    print("ancestries are: ")
    print(ancestries)

    out1 <- algorithm_1snp_wrapper(i = 1,
                                  alleles = alleles,
                                  ancestries = ancestries)

    print("algorithm_1snp_wrapper(i = 1,
                                  alleles = alleles,
                                  ancestries = ancestries)")
    print(out1)

    out2 <- algorithm_1snp_wrapper(i = 2,
                                  alleles = alleles,
                                  ancestries = ancestries)

    print("algorithm_1snp_wrapper(i = 2,
                                  alleles = alleles,
                                  ancestries = ancestries)")
    print(out2)
})
