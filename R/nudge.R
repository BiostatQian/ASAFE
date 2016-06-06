# nudge: Nudge a probability estimate away from 0
# if it is on the boundary (either 0 or 1).
#
# INPUTS:
# 1) Numeric value p: A probability estimate
# 2) Epsilon: Amount by which to nudge p away from 0 or 1
#
# OUTPUTS:
# 1) Nudged p (i.e. p + epsilon if input p = 0, p - epsilon, if input p = 1)

nudge <- function(p, epsilon = 10^-8){
    # print("nudge")

    if(p == 0){
        #print("p = 0")
        out <- p + epsilon
    }else if(p == 1){
        #print("p = 1")
        out <- p - epsilon
    }else{
        #print("p not 0 or 1")
        out <- p
    }

    return(out)
}
