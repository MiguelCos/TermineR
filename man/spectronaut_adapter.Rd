% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adapters.R
\name{spectronaut_adapter}
\alias{spectronaut_adapter}
\title{spectronaut_adapter}
\format{
A data frame with at least 4 columns
\describe{
\item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
\item{nterm_modif}{N-terminal modification. Only Acetyl annotation is available in the current version}
\item{peptide}{Peptide sequence}
\item{protein}{Protein ID based on Uniprot Accession}
}
}
\usage{
spectronaut_adapter(path_to_file, proteotypic = TRUE)
}
\arguments{
\item{path_to_file}{path to file with Spectronaut search results in tsv format}

\item{proteotypic}{keep only proteotypic peptides. Default is TRUE.}
}
\description{
Process report file from Spectronaut search results for further analysis with TermineR
}
\details{
Process report file from Spectronaut search results for further analysis with TermineR
}
\author{
Miguel Cosenza-Contreras
}
