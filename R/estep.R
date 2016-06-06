# E-Step
#
# INPUTS:
#
# 1) Number of people in each of 18 observed categories
# n = [n1, ..., n18].
# Expected class: Numeric.
#
# 2) Latest estimate of p = [p_GA, p_GE, p_GN, p_TA, p_TE, p_TN],
# i.e. p = [p_00, p_01, p_02, p_10, p_11, p_12],
# where we consider 0/1 = G/T and 0/1/2 = A/E/N.
# Indices of p are [1, 2, 3, 4, 5, 6].
# Expected class: Numeric.
#
# OUTPUTS:
#
# 1) Expectations of number of people in each of 21 complete
# (g1,a1)/(g2,a1) categories, m = 1, ..., 21.
# Expected class: Numeric.

estep <- function(mprime, p){
    # print("EStep")

    # Output vector of expectations
    m <- numeric(21)

    # 5/25/2016: 
    # Implementing the supplement directly,
    # as given in column 4 of "Table 2" of the supplement.

    m[1] <- mprime[1]
    m[2] <- mprime[2]
    m[3] <- mprime[3]
    m[4] <- mprime[4]
    m[5] <- mprime[5]
    m[6] <- mprime[6]
    m[7] <- mprime[7]
    m[8] <- mprime[8] * (p[1] * p[5]) / (p[1] * p[5] + p[2] * p[4])
    m[9] <- mprime[8] * (p[2] * p[4]) / (p[1] * p[5] + p[2] * p[4])
    m[10] <- mprime[9] * (p[1] * p[6]) / (p[1] * p[6] + p[3] * p[4])
    m[11] <- mprime[9] * (p[3] * p[4]) / (p[1] * p[6] + p[3] * p[4])
    m[12] <- mprime[10]
    m[13] <- mprime[11] * (p[2] * p[6]) / (p[2] * p[6] + p[3] * p[5])
    m[14] <- mprime[11] * (p[3] * p[5]) / (p[2] * p[6] + p[3] * p[5])
    m[15] <- mprime[12]
    m[16] <- mprime[13]
    m[17] <- mprime[14]
    m[18] <- mprime[15]
    m[19] <- mprime[16]
    m[20] <- mprime[17]
    m[21] <- mprime[18]

    return(m)
}
