\name{algorithm_1snp_wrapper}
\alias{algorithm_1snp_wrapper}
\title{
Wrapper for function algorithm_1snp
}
\description{
Applies the function \link{algorithm_1snp} to a particular bi-allelic marker's data stored in matrices alleles and ancestries.
Can be used to apply \link{algorithm_1snp} to a set of bi-allelic markers.
}
\usage{
algorithm_1snp_wrapper(i, alleles, ancestries)
}
\arguments{
  \item{i}{Index of marker in matrices alleles and ancestries.
  This is the index for the marker that we want to apply function \link{algorithm_1snp} to.
}
  \item{alleles}{
  Rows: Alleles for individuals' chromosomes ordered
e.g. ADM1, ADM1, ..., ADM250, ADM250, where ADMi is the ID for the i-th individual.
Cols: Bi-allelic markers.
Each allele is either 0 or 1.
This is a numeric matrix.
}
  \item{ancestries}{
  Rows: Ancestries for individuals' chromosomes ordered
e.g. ADM1, ADM1, ..., ADM250, ADM250, where ADMi is the ID for the i-th individual.
Cols: Bi-allelic markers.
Each ancestry is either 0, 1, or 2.
This is a numeric matrix.
}
}
\details{Markers in matrix alleles should be in 1-to-1 correspondence
with markers in matrix ancestries. Markers in both matrices should
be in the same order.
}
\value{
A character vector with 4 elements.
The first element is the Marker ID
of the i-th marker in matrices alleles and ancestries.
The next 3 elements are
ancestry-specific allele frequency estimates of
P(Allele 1 | Ancestry 0), P(Allele 1 | Ancestry 1), and P(Allele 1 | Ancestry 2),
for the i-th marker in matrices alleles and ancestries.
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

adm_ancestries_test <- head(adm_ancestries)

# adm_genotypes_test is a matrix with
# Rows: Markers
# Columns: Marker ID, individuals' genotypes (a1/a2)
# (e.g. ADM1, ADM2, ADM3, and etc.)

adm_genotypes_test <- head(adm_genotypes)

# Making the rsID column row names
row.names(adm_ancestries_test) <- adm_ancestries_test[,1]
row.names(adm_genotypes_test) <- adm_genotypes_test[,1]

adm_ancestries_test <- adm_ancestries_test[,-1]
adm_genotypes_test <- adm_genotypes_test[,-1]

# alleles_list is a list of lists.
# Outer list elements correspond to bi-allelic markers.
# Inner list elements correspond to 250 people's alleles with no delimiter separating alleles.
alleles_list <- apply(X = adm_genotypes_test, MARGIN = 1, FUN = strsplit, split = "/")

# Creates a matrix: Number of alleles (ADM1, ADM1, ..., ADM250, ADM250) x (bi-allelic markers)
alleles_unlisted <- sapply(alleles_list, unlist)

# Change elements of the matrix to numeric, producing a matrix:
# Number of alleles (ADM1, ADM1, ..., ADM250, ADM250) x (bi-allelic markers).
alleles <- apply(X = alleles_unlisted, MARGIN = 2, as.numeric)

# Apply EM algorithm to first bi-allelic marker
algorithm_1snp_wrapper(i = 1, alleles, adm_ancestries_test)

# Apply the EM algorithm to each bi-allelic marker to obtain
# ancestry-specific allele frequency estimates for all markers in
# matrices alleles and ancestries.

adm_estimates_test <- sapply(X = 1:ncol(alleles), FUN = algorithm_1snp_wrapper,
                        alleles = alleles, ancestries = adm_ancestries_test)

adm_estimates_test

}

