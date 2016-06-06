\name{algorithm_1snp}
\alias{algorithm_1snp}
\title{
Estimate ancestry-specific allele frequencies for 1 marker (e.g. a SNP) from individuals' alleles and ancestries at this marker.
}
\description{
Take in genotypes (possibly unphased with respect to each other) and ancestries (possibly unphased with respect to each other) for all individuals at 1 marker to create the marker's vector of observed data category counts, and then call the function em() on that vector of counts, to obtain ancestry-specific allele frequency estimates for that marker.
}
\usage{
algorithm_1snp(alleles_1, ancestries_1)
}
\arguments{
  \item{alleles_1}{Vector of alleles for each individual's 2 chromosomes, with chromosomes for the same individual consecutive. Each allele is either 0 or 1. This is a numeric vector.

Example: If there are 250 admixed individuals, the alleles might be ordered like so: ADM1, ADM1, ADM2, ADM2, ..., ADM250, ADM250, where ADMi is the ID for the i-th individual.
}
  \item{ancestries_1}{Vector of ancestries for each individual's 2 chromosomes,
with chromosomes for the same individual consecutive. Each ancestry is either 0, 1, or 2. This is a numeric vector.

Example: If there are 250 admixed individuals, the ancestries might be ordered like so: ADM1, ADM1, ADM2, ADM2, ..., ADM250, ADM250, where ADMi is the ID for the i-th individual.
}
}
\value{
Ancestry-specific allele frequency estimates of [P(Allele 1| Ancestry 0), P(Allele 1 | Ancestry 1), P(Allele 1 | Ancestry 2)] from the EM Algorithm. This a numeric vector with 3 entries.
}
\author{Qian Zhang
}
\examples{

adm_ancestries_test <- head(adm_ancestries)
adm_genotypes_test <- head(adm_genotypes)

# adm_ancestries_test is a matrix with
# Rows: Markers
# Columns: Marker ID, individuals' chromosomes' ancestries
# (e.g. ADM1, ADM1, ADM2, ADM2, and etc.)

# adm_genotypes_test is a matrix with
# Rows: Markers
# Columns: Marker ID, individuals' genotypes (a1/a2)
# (e.g. ADM1, ADM2, ADM3, and etc.)

# Make the rsID column row names
row.names(adm_ancestries_test) <- adm_ancestries_test[,1]
row.names(adm_genotypes_test) <- adm_genotypes_test[,1]

adm_ancestries_test <- adm_ancestries_test[,-1]
adm_genotypes_test <- adm_genotypes_test[,-1]

# alleles_list is a list of lists.
# Outer list elements correspond to SNPs.
# Inner list elements correspond to 250 individuals's alleles with no delimiter separating alleles.

alleles_list <- apply(X = adm_genotypes_test, MARGIN = 1,
                        FUN = strsplit, split = "/")

# Creates a matrix: Number of alleles
# (ADM1, ADM1, ..., ADM250, ADM250) x (SNPs)

alleles_unlisted <- sapply(alleles_list, unlist)

# Change elements of the matrix to numeric, producing a matrix:
# Number of alleles (ADM1, ADM1, ..., ADM250, ADM250) x (SNPs).

alleles <- apply(X = alleles_unlisted, MARGIN = 2, as.numeric)

# Perform the EM algorithm on the first SNP in the data, obtaining estimates for
# P(Allele 1 | Ancestry 0), P(Allele 1 | Ancestry 1), P(Allele 1 | Ancestry 2)

estimates <- algorithm_1snp(alleles[,1], adm_ancestries_test[1,])

estimates
}

