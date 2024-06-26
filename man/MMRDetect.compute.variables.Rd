% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMRDetect.compute.variables.R
\name{MMRDetect.compute.variables}
\alias{MMRDetect.compute.variables}
\title{Compute variables for MMRDetect}
\usage{
MMRDetect.compute.variables(
  sub_cat,
  indel_cat,
  tissue_type,
  MMR_subsig96 = MMRKO_subsig,
  MMR_sig_indel = MMRKO_indelsig,
  tissue_subsig96 = PancanSig
)
}
\arguments{
\item{sub_cat}{substitution catalogue}

\item{indel_cat}{indel catalogue}

\item{tissue_type}{tissue type}

\item{MMR_subsig96}{mismatch repair gene knockout substitution signatures}

\item{MMR_sig_indel}{mismatch repair gene knockout indel signatures}

\item{tissue_subsig96}{tissue-specific substitution signatures}
}
\value{
variables
}
\description{
Compute variables for MMRDetect
}
