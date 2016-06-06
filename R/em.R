##### EM
#
# INPUTS:
#
# 1) n = [n1, ..., n18] of observed data categories counts,
# i.e. number of people in each of 18 categories.
#
# 2) epsilon: 
# Difference between estimates such that we consider the algorithm converged.
# Default value: 10^-8.
#
# 3) iteration_cap: Max number of iterations we'll allow EM to go for.
# Default value: 1000.
#
# OUTPUTS:
#
# 1) Vector of conditional ancestry-specific allele frequencies
# [P(Alelle 1 | Ancestry 0), P(Allele 1 | Ancestry 1), 
# P(Allele 1 | Ancestry 2)].
#
# DETAILS:
#
# 0) Start with initial guess for p_{k=0} (arbitrary). k indexes the iteration.
# k = 1, 1) Do iteration k = 1's E-Step. Start with p_{k = 0}. Denote pk.
# k = 1, 2) Do iteration k = 1's M-Step. End with p_{k = 1}. Denote pk1.
# k = 2, 1) Do iteration k = 2's E-Step. Start with p_{k = 1}. Denote pk.
# k = 2, 2) Do iteration k = 2's M-Step. End with p_{k = 2}. Denote pk1.
# Etc.
# p_{k} is the estimate that you end iteration k with.
# Stop iterating when ||pk - pk1|| < epsilon,
# i.e. when the last 2 estimates' are close to each other.

em <- function(n, epsilon = 10^-8, iteration_cap = 1000){
    # message("EM")

    # Initialize iteration index k = 0.
    k <- 0

    # Initialize vector of EM ancestry-specific allele 0/1 frequencies
    # for ancestries 0, 1, and 2. pk = [p00, p01, p02, p10, p11, p12].
    pk <- rep(1/6, times = 6) # Last estimate (e.g. iteration k = 0)
    # message(k)

    # k = 1, 1) Iteration k = 1's E-Step
    m_k <- estep(mprime = n, p = pk)

    # k = 1, 2) Iteration k = 1's M-Step
    pk1 <- mstep(m = m_k) # New estimate (e.g. iteration k = 1)

    # Update k to say which iteration just finished
    k <- 1
    #message(k)

    # Continue iterating the algorithm so long as difference
    # between estimate at beginning and at end
    # of latest iteration are very different, i.e. ||pk+1 - pk|| > epsilon,
    # and we have not reached the cap iteration_cap 
    # on the number of iterations.

    while(sqrt(sum((pk1 - pk)^2)) > epsilon){

        # Update last estimate pk.
        # This is the pk estimate at the beginning of this iteration k,
        # i.e. going into this iteration k.

        pk <- pk1
        m_k <- estep(mprime = n, p = pk)

        # This is the new pk estimate at the end of this iteration k,
        # which will be used as the estimate at the beginning of 
        # iteration k + 1

        pk1 <- mstep(m = m_k)

        # Update k to reflect the iteration of the algorithm that just finished
        k <- (k+1)
        # message(k)

        # If we have reached the cap on number of iterations, print warning
        if(k == iteration_cap){
            message(paste("Iteration cap was ", iteration_cap))
            warning('Max number of allowed iterations reached. 
                     Outputting latest result.')
            break # Exit while loop
        }
    }

    # pk1 should now give EM estimates for p = [p00, p01, p02, p10, p11, p12].
    # message(paste("Allele 1 frequencies conditional on 
    # ancestries 0, 1, 2 found in ", k, iterations"))

    # pk1_conditional = 
    # [P(allele 1 | anc 0), P(allele 1 | anc 1), P(allele 1 | anc 2)].

    pk1_conditional <- rep(0, times = 3)

    pk1_conditional[1] <- pk1[4] / (pk1[1] + pk1[4])
    pk1_conditional[2] <- pk1[5] / (pk1[2] + pk1[5])
    pk1_conditional[3] <- pk1[6] / (pk1[3] + pk1[6])

    return(pk1_conditional)
}
