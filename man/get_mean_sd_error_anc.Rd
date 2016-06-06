\name{get_mean_sd_error_anc}
\alias{get_mean_sd_error_anc}
\title{
Get summary statistics on errors
for markers with a certain ancestry's
true allele 1 frequency in an allele frequency bin.
}
\description{
For each (min_bin, max_bin] allele frequency bin
i.e. (0,0.2], (0.2, 0.4], (0.4, 0.6], (0.6, 0.8], and (0.8, 1.0],
get mean and SD of errors
(where error = estimated allele 1 frequency - true allele 1 frequency)
for markers with a certain ancestry's true allele 1 frequency in the bin.
}
\usage{
get_mean_sd_error_anc(ancestry, estimates, truth)
}
\arguments{
  \item{ancestry}{Ancestry i.e. ancestral population, coded 0, 1, or 2}
  \item{estimates}{Estimated allele frequencies

Rows: Markers

Columns: 1) Marker ID,
2-4) Ancestry-specific allele 1 frequency estimates for ancestries 0, 1, and 2
}
  \item{truth}{True allele frequencies.

Rows: Markers

Columns: 1) Marker ID,
2-4) Ancestry-specific allele 1 frequency true values for ancestries 0, 1, and 2
}
}
\value{
  \item{out}{A 2 x 5 matrix with

Rows: (1) mean of errors, (2) sd of errors

Cols: True ancestry-specific allele 1 frequency bins
(0,0.2], (0.2, 0.4], (0.4, 0.6], (0.6, 0.8],
(0.8,1.0]

}
}

\author{Qian Zhang
}
\examples{

# Simulated ancestries for 56,003 markers
# Rows: 56,003 markers
# Columns: 1) Marker ID, 2-4) Ancestry 0, 1, 2 allele 1 frequency estimates
dim(adm_estimates)

# ancestral_freqs
# Rows: 56,003 Markers
# Columns: 1) Marker ID,
# 2-4) Ancestry 0, 1, 2 allele 1 frequency true values
dim(ancestral_freqs)

# Mean and SD of errors for markers with
# true ancestry 0 allele 1 frequency in bins
# (0, 0.2], (0.2, 0.4], (0.4, 0.6], (0.6, 0.8], (0.8, 1.0].

get_mean_sd_error_anc(ancestry = 0, estimates = adm_estimates, truth = ancestral_freqs)
}

