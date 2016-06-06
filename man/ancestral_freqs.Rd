\name{ancestral_freqs}
\docType{data}
\alias{ancestral_freqs}
\title{True ancestry-specific allele frequencies for 56,003 markers}
\description{
True ancestry-specific allele frequencies computed
from simulated phased genotypes
and simulated phased ancestries
}
\usage{ancestral_freqs}
\format{

A 56,003 x 4 matrix with the following rows, columns, and entries:

    \enumerate{
    \item Rows: 1 row per bi-allelic marker,
    with alleles arbitrarily labeled 0 and 1
    \item Columns: Marker ID, followed by ancestries 0, 1, and then 2
    \item Entries: The (i,j)-th entry
    when j is a colum corresponding to an ancestry,
    gives the true frequency of allele 1 amongst admixed chromosomes
    of that ancestral origin at the i-th marker.
    }
}
\source{Simulated genetic data}
\references{
Zhang QS, Browning BL, and Browning SR (2016) ''Ancestry Specific Allele Frequency Estimation.'' Bioinformatics.
}
\author{Qian Zhang}

