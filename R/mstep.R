# MStep: Perform M-step of EM algorithm
#
# INPUTS:
# (1) Expected complete data category counts m = [m1, ..., m21].
# (2) epsilon: Amount by which you push a probability estimate
# above 0 or below 1. Default value is 10^-8.
#
# OUTPUTS:
# (1) Updated estimates of probabilities p =
# [p00, p01, p02, p10, p11, p12],
# i.e. [p_GA, p_GE, p_GN, p_TA, p_TE, p_TN].
# i.e. [p1, p2, p3, p4, p5, p6].
#
# DETAILS: 
# If a probability has an estimate that is identically 0
# (on a boundary, since a probabilty can only be in [0,1]),
# push the estimate off the boundary,
# so that you don't get divide by 0 error in the EStep.
# I'll also nudge the estimate below 1,
# because combined with nudging an estimate above 0,
# this will help to ensure that the sum of frequency estimates is 1.
# 
# Indexing of complete data category vector m 
# follows ASAFE Supplement pages 13-14.
# Code's [p[1], ..., p[6]] corresponds to 
# Supplement's [p_ga, g in {0,1}, a in {A,E,N}] page 15. 

mstep <- function(m){
    # print("MStep")

    # number of observations (i.e. individuals) n = Sum of m vector
    n <- sum(m)
    denominator <- 2*n

    # Define output vector
    p <- numeric(6)

    # Numerator for p_GA i.e. p_00 i.e. p[1].
    # In terms of supplement indices:
    # N0A = 2*m[1] + m[2] + m[3] + m[7] + m[8] + m[10].

    p[1] <- 2*m[1] + m[2] + m[3] + m[7] + m[8] + m[10]
    p[1] <- p[1] / denominator

    # Numerator for p_GE i.e. p_01 i.e. p[2].
    # In terms of supplement indices:
    # N0E = m[2] + 2*m[4] + m[5] + m[9] + m[12] + m[13].

    p[2] <- m[2] + 2*m[4] + m[5] + m[9] + m[12] + m[13]
    p[2] <- p[2] / denominator

    # Numerator for p_GN i.e. p_02 i.e. p[3].
    # In terms of supplement indices:
    # N0N = m[3] + m[5] + 2*m[6] + m[11] + m[14] + m[15].

    p[3] <- m[3] + m[5] + 2*m[6] + m[11] + m[14] + m[15]
    p[3] <- p[3] / denominator

    # Numerator for p_TA i.e. p_10 i.e. p[4].
    # In terms of supplement indices:
    # N1A = m[7] + m[9] + m[11] + 2*m[16] + m[17] + m[18].

    p[4] <- m[7] + m[9] + m[11] + 2*m[16] + m[17] + m[18]
    p[4] <- p[4] / denominator

    # Numerator for p_TE i.e. p_11 i.e. p[5].
    # In terms of supplement indices:
    # N1E = m[8] + m[12] + m[14] + m[17] + 2*m[19] + m[20].

    p[5] <- m[8] + m[12] + m[14] + m[17] + 2*m[19] + m[20]
    p[5] <- p[5] / denominator

    # Numerator for p_TN i.e. p_12 i.e. p[6].
    # In terms of supplement indices:
    # N1N = m[10] + m[13] + m[15] + m[18] + m[20] + 2*m[21].

    p[6] <- m[10] + m[13] + m[15] + m[18] + m[20] + 2*m[21]
    p[6] <- p[6] / denominator

    # Nudge each probability epsilon above 0 if it is 0,
    # or epsilon below 1 if it is 1.
    # Epsilon is built into the nudge function.

    out <- sapply(p, nudge)

    return(out)
}
