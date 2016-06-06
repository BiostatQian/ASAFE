\name{adm_estimates}
\docType{data}
\alias{adm_estimates}
\title{Ancestry-specific allele frequency estimates for 56,003 markers}
\description{
Ancestry-specific allele frequency estimates computed
by applying ASAFE's EM algorithm
to simulated unphased genotypes
(ignoring phase information from simulated phased genotypes)
and simulated phased ancestries.
}
\usage{adm_estimates}
\format{

A 56,003 x 4 matrix with the following rows, columns, and entries:

    \enumerate{
    \item Rows: 1 row per bi-allelic marker,
    with alleles arbitrarily labeled 0 and 1
    \item Columns: Marker ID, followed by ancestries 0, 1, and then 2
    \item Entries: The (i,j)-th entry when j is a column
    corresponding to an ancestry,
    gives the estimated frequency of allele 1 amongst admixed chromosomes
    of that ancestral origin at the i-th marker.
    }
}
\source{This is the result of applying ASAFE functions to
simulated genotype and ancestry data. See ASAFE vignette for details.}
\references{
Zhang QS, Browning BL, and Browning SR (2016) ''Ancestry Specific Allele Frequency Estimation.'' Bioinformatics.
}
\author{Qian Zhang}


