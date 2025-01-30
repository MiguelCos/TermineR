#' @title Mouse protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'
#' @format A data frame
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

#' @title Human protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'
#' @format A data frame
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
"human_uniprot_processing"

#' @title Human protein processing annotation from Uniprot API. Swissprot including isoforms.
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'
#' @format A data frame
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
#' @source Uniprot API (queried on August 2024)
"human_and_isoforms_uniprot_processing"

#' @title Mendicato trucantula protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
"mendicato_trucantula_uniprot_processing"

#' @title Rhizobium melitoli protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
"rhizobium_melitoli_uniprot_processing"

#' @title Sus scrofa protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
#' @source Uniprot API (queried on July 2024)
"pig_uniprot_processing"

#' @title UniMOd ID to modification name mapping 
#'
#' @description Tabular annotation of the numerical ID from UniMod and their associated modifications
#'   
#' @format A data frame
#' \describe{
#'  \item{id}{UniMod ID}
#'  \item{id_nr}{Numerical UniMod ID}
#'  \item{name}{Descriptive name of the modification}
#' }
#' @source Uniprot API (queried on August 2024)
"unimod_id_to_name_mapping"

#' @title Arabidopsis thaliana protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
#' @source Uniprot API (queried on September 2024)
"arabidopsis_uniprot_processing"

#' @title Escherichia coli protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
#' @source Uniprot API (queried on September 2024)
"ecoli_uniprot_processing"


#' @title Rattus norvegicus protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
#' @source Uniprot API (queried on September 2024)
"rat_uniprot_processing"

#' @title Saccharomyces cerevisiae protein processing annotation from Uniprot API
#'
#' @description Tabular annotation of processing information as obtained from Uniprot API
#'   
#' @format A data frame
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
#' @source Uniprot API (queried on September 2024)
"yeast_uniprot_processing"