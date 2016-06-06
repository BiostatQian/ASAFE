\name{adm_genotypes}
\docType{data}
\alias{adm_genotypes}
\title{Genotypes of 250 admixed individuals at 56,003 markers}
\description{

This matrix adm_genotypes stores the full data set of simulated phased genotypes (meaning genotypes at different markers are phased with respect to each other) that was used in the ASAFE paper. adm_genotypes contains genotypes at 56,003 markers for 250 admixed individuals.

For each individual at a marker, the genotype is also phased with respect to the ancestry pair given in \link{adm_ancestries}, so that true ancestry-specific allele frequencies can be calculated from adm_genotypes and \link{adm_ancestries} to obtain \link{ancestral_freqs}, by overlaying ancestries on genotypes. The ASAFE EM algorithm does not assume that ancestries and genotypes at the same marker are phased with respect to each other, or that ancestries at different markers are phased with respect to each other, or that genotypes at different markers are phased with respect to each other, and provides estimates of true ancestry-specific allele frequencies.

}
\usage{adm_genotypes}
\format{

    A 56,003 x 251 matrix with the following rows, columns, and entries:

    \enumerate{
    \item Rows: 1 row per bi-allelic marker,
    with alleles arbitrarily labeled 0 and 1
    \item Columns: First column is Marker ID. Following columns consist of 1 column per individual. Individuals should be listed in the same order in the genotype matrix adm_genotypes as in the ancestry matrix adm_ancestries.
    \item Entries: For an entry that is not in the Marker ID column,
an entry can take value 0/0, 0/1, 1/0, or 1/1, where 0 and 1 are arbitrary labels for a bi-allelic SNP's two alleles. A slash ''/'' indicates an unphased genotype, so 0/1 and 1/0 are the same unphased genotype. It just so happens that this particular adm_genotypes matrix contains phased genotypes, despite the presence of slashes. However, the ASAFE algorithm does not assume that genotypes are phased with respect to each other.
}

}
\source{Simulated genetic data}
\references{
Zhang QS, Browning BL, and Browning SR (2016) ''Ancestry Specific Allele Frequency Estimation.'' Bioinformatics.
}
\author{Qian Zhang}



