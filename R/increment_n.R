# increment_n: Contribute an individual's unphased genotype,
# phased ancestry data to the observed data category counts.
#
# INPUTS:
#
# 1) Vector of observed data category counts n = [n1, n2, ..., n18].
#
# 2) 2-vector of the person's ancestries (a1, a2)
# where a1, a2 are each either 0, 1, or 2.
#
# 3) 2-vector of the person's unphased genotype alleles (g1, g2),
# where g1, g2 are either 0 or 1. (g1, g2) suggests order,
# but in these comments it means the same thing as g1/g2
# just with the slash / removed.
#
# OUTPUTS:
# 1) Updated vector of observed data category counts n = [n1, n2, ..., n18].
#
# DETAILS: 
# In the ASAFE paper supplement, vector n is called m_prime.

increment_n <- function(n, anc_ind, geno_ind){
    # print("increment_n")

    anc_copy1 <- anc_ind[1]
    anc_copy2 <- anc_ind[2]

    allele1_copies <- sum(geno_ind)

    # alleles1_copies is the number of copies of the type 1 allele.
    # Ancestries 0, 1, and 2 are called A, E, and N in the supplement.
    # Indexing of vector n follows the observed data category definitions
    # given in the 4th column from the left 
    # of Supplementary Table 1 of the ASAFE paper.

    if(allele1_copies == 0){
        if(anc_copy1 == 0 & anc_copy2 == 0){
            n[1] <- n[1] + 1
        }else if(anc_copy1 == 1 & anc_copy2 == 1){
            n[4] <- n[4] + 1
        }else if(anc_copy1 == 2 & anc_copy2 == 2){
            n[6] <- n[6] + 1
        }else if((anc_copy1 == 1  & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 1)){
            n[2] <- n[2] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 2)){
            n[3] <- n[3] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 1) || (anc_copy1 == 1 & anc_copy2 == 2)){
            n[5] <- n[5] + 1
        }else{
            print("WARNING: Impossible ancestry combo")

            print("Number of copies of allele 1: ")
            print(allele1_copies)

            print("Ancestries are: ")
            print(c(anc_copy1, anc_copy2))

            print("Quitting program ")

            # Quit program
            quit()
        }
    }else if(allele1_copies == 1){
        if(anc_copy1 == 0 & anc_copy2 == 0){
            n[7] <- n[7] + 1
        }else if(anc_copy1 == 1 & anc_copy2 == 1){
            n[10] <- n[10] + 1
        }else if(anc_copy1 == 2 & anc_copy2 == 2){
            n[12] <- n[12] + 1
        }else if((anc_copy1 == 1  & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 1)){
            n[8] <- n[8] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 2)){
            n[9] <- n[9] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 1) || (anc_copy1 == 1 & anc_copy2 == 2)){
            n[11] <- n[11] + 1
        }else{
            print("WARNING: Impossible ancestry combo")

            print("Number of copies of allele 1: ")
            print(allele1_copies)

            print("Ancestries are: ")
            print(c(anc_copy1, anc_copy2))

            print("Quitting program ")

            # Quit program
            quit()
        }
    }else if(allele1_copies == 2){
        if(anc_copy1 == 0 & anc_copy2 == 0){
            n[13] <- n[13] + 1
        }else if(anc_copy1 == 1 & anc_copy2 == 1){
            n[16] <- n[16] + 1
        }else if(anc_copy1 == 2 & anc_copy2 == 2){
            n[18] <- n[18] + 1
        }else if((anc_copy1 == 1  & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 1)){
            n[14] <- n[14] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 0) || (anc_copy1 == 0 & anc_copy2 == 2)){
            n[15] <- n[15] + 1
        }else if((anc_copy1 == 2 & anc_copy2 == 1) || (anc_copy1 == 1 & anc_copy2 == 2)){
            n[17] <- n[17] + 1
        }else{
            print("WARNING: Impossible ancestry combo")

            print("Number of copies of allele 1: ")
            print(allele1_copies)

            print("Ancestries are: ")
            print(c(anc_copy1, anc_copy2))

            print("Quitting program ")

            # Quit program
            quit()
        }
    }else{
        print("WARNING: Impossible number of copies of allele 1")
    }

    return(n)
}

