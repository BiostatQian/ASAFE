# INPUTS:
# (1) p0: Ancestry 0 (i.e. ancestral population 0) allele 1 frequency 
# (2) p1: Ancestry 1 (i.e. ancestral population 1) alelle 1 frequency
# (3) p2: Ancestry 2 (i.e. ancestral population 2) allele 1 frequency
#
# (4) ancestries_matrix_true: True ancestries
# Rows: Individuals' alleles (1, ..., 2*n, when there are n individuals) 
# Cols: Bi-allelic markers (e.g. SNPs)
#
# (5) ancestries_matrix_estimated: Ancestries that ASAFE uses for 
# estimation of ancestry specific allele frequencies. 
# May or may not be true ancestries.
# Rows: Individuals' alleles (1, ..., 2*n, when there are n individuals) 
# Cols: Bi-allelic markers (e.g. SNPs) 
#
# WHAT THE FUNCTION DOES:
# 
# (1) Given ancestral populations' allele 1 frequencies p0, p1, and p2,
# draw alleles given true ancestries, 
# (2) Run ASAFE with drawn alleles and ancestries in 
# ancestries_matrix_estimated to obtain 
# P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2) estimates. 
# (3) Calculate true 
# P(Allele 1 | Anc 0), P(Allele 1 | Anc 1), P(Allele 1 | Anc 2)
# from from alleles and true ancestries in ancestries_matrix_true.
# (4) Calculate error in estimating frequencies (estimates - true values). 
# 
# OUTPUTS:
# (1) (Number of ancestries i.e. 3) x (Number of bi-allelic markers) 
# matrix of errors. Each entry gives error in 
# P(Allele 1 | Ancestry a) as estimate - true value.
 
get_errors_1_scenario <- function(p0, p1, p2, ancestries_matrix_true,
                                                ancestries_matrix_estimated){

    # Alleles matrix.
    # - Rows: Individuals' alleles 
    # (1, ..., 2*n, when n = number of individuals). 
    # - Cols: Markers.

    alleles_matrix <- matrix(mapply(FUN = function(anc){
                                        draw_allele_given_anc(anc, 
                                                              freqs = c(p0, p1, p2))
                                    }, ancestries_matrix_true),
                                    nrow = nrow(ancestries_matrix_true),
                                    ncol = ncol(ancestries_matrix_true))

    # Outputs a (Number of ancestries) x 
    # (Number of markers) matrix of 
    # estimated ancestry-specific allele frequencies per marker.
    # Row 1: P(Allele 1 | Anc 0) estimate.
    # Row 2: P(Allele 1 | Anc 1) estimate.
    # Row 3: P(Allele 1 | Anc 2) estimate.
    # Cols: Markers.

    estimates <- sapply(X = 1:ncol(alleles_matrix),
                        FUN = function(marker_col){
                            algorithm_1snp(alleles_1 = alleles_matrix[,marker_col],
                                           ancestries_1 = ancestries_matrix_estimated[,marker_col])
                        })

    # Outputs a (Number of ancestries) x 
    # (Number of markers) matrix of 
    # true ancestry-specific allele frequencies per marker.
    # Row 1: P(Allele 1 | Anc 0) true value.
    # Row 2: P(Allele 1 | Anc 1) true value.
    # Row 3: P(Allele 1 | Anc 2) true value.
    # Cols: Markers.
 
    true_freqs <- sapply(X = 1:ncol(alleles_matrix),
                         FUN = function(marker_col){
                             get_true_freqs_1snp(alleles_1 = alleles_matrix[,marker_col],
                                                 ancestries_1 = ancestries_matrix_true[,marker_col])
                         })

    # Outputs a (Number of ancestries) x 
    # (Number of markers) matrix of errors
    # (estimated allele freq - true allele freq),
    # for each marker.
   
    errors <- (estimates - true_freqs)


    out <- errors

    return(out)

}
