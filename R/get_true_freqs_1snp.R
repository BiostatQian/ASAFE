# Get true ancestry-specific allele frequencies

get_true_freqs_1snp <- function(alleles_1, ancestries_1){

    indices_anc_0 <- which(ancestries_1 == 0)
    anc0_freq <- sum(alleles_1[indices_anc_0]) / length(indices_anc_0)

    indices_anc_1 <- which(ancestries_1 == 1)
    anc1_freq <- sum(alleles_1[indices_anc_1]) / length(indices_anc_1)

    indices_anc_2 <- which(ancestries_1 == 2)
    anc2_freq <- sum(alleles_1[indices_anc_2]) / length(indices_anc_2)

    return(c(anc0_freq, anc1_freq, anc2_freq))
}
