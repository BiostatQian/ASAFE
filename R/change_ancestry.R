# Change ancestry. 
# Decide to change ancestry with probability error_rate.
# If decide to change ancestry, 
# then sample the other two ancestries with 50/50 chance.

change_ancestry <- function(anc, error_rate){

    ind_change <- sample(x = 0:1, size = 1, 
                         prob = c(1 - error_rate, error_rate))

    if(ind_change == 1){
        message("ind_change == 1")
        out <- sample_ancestry(anc)
    }else{
        out <- anc
    }

    return(out)
}

