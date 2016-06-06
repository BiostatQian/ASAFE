# INPUTS
#
# 1) i: Index of SNP to analyze
#
# 2) alleles:
# - Rows: Alleles ordered e.g. ADM1, ADM1, ..., ADM250, ADM250.
# - Cols: Names of markers (e.g. SNPs)
# - Matrix is assumed to have column names that are marker ID names.
#
# 3) ancestries:
# - Rows: Ancestries ordered e.g. ADM1, ADM1, ..., ADM250, ADM250.
# - Cols: Names of markers (e.g. SNPs)
#
# OUTPUTS
#
# 1) Ancestry-specific allele frequencies for this SNP:
# P(1|A), P(1|E), P(1|N).

algorithm_1snp_wrapper <- function(i, alleles, ancestries){
    rsid <- colnames(alleles)[i]
    estimates <- algorithm_1snp(alleles[,i], ancestries[i,])
    out <- c(rsid, estimates)

    return(out)
}
