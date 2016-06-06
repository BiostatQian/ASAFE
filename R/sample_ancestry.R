###### Add in ancestry-calling error.
# Change each ancestry to a different ancestry with equal probability
# of being one of the other two ancestries.

sample_ancestry <- function(anc){

    if(anc == 0){
        out <- sample(x = c(1,2), size = 1, prob = c(0.5, 0.5))
    }else if(anc == 1){
        out <- sample(x = c(0,2), size = 1, prob = c(0.5, 0.5))
    }else if(anc == 2){
        out <- sample(x = c(0,1), size = 1, prob = c(0.5, 0.5))
    }else{
        stop("anc should be 0, 1, or 2")
    }

    return(out)
}
