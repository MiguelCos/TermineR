#' @title Protein processing annotation from Uniprot API
#'
#' @description Annotation table of samples and conditions for use case experiment of
#'      polycystic kidney disease in mice.
#'
#' @format A data frame with 9 columns anmd 29256 rows
#' \describe{
#'  \item{type}{Type of processing feature. One of SIGNAL, CHAIN, INIT_MET, TRANSIT, PEPTIDE}
#'  \item{description}{Description of processing feature}
#'  \item{begin}{begin of the processing feature in protein sequence}
#'  \item{end}{end of the processing feature in protein sequence}
#'  \item{length}{length of the processing feature within protein sequence}
#'  \item{accession}{Uniprot accession ID associated with the feature}
#'  \item{entryName}{Uniprot Entry name ID associated with the feature}
#'  \item{taxid}{Organism tax id}
#'  \item{order}{order}
#' }
#' @source Uniprot API (queried on April 2024)
"mouse_uniprot_processing"
