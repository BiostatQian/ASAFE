# Draw an allele given ancestry 
# and vector of ancestry-specific allele frequencies.

draw_allele_given_anc <- function(anc, freqs){

    if(anc == 0){
        allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[1], freqs[1]))
    }else if(anc == 1){
        allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[2], freqs[2]))
    }else if(anc == 2){
        allele <- sample(x = 0:1, size = 1, prob = c(1 - freqs[3], freqs[3]))
    }else{
        stop('anc must be 0, 1, or 2')
    }

    return(allele)
}

