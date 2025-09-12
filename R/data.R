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

#' @title Caenorhabditis elegans protein processing annotation from Uniprot API
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
#' @source Uniprot API (queried on December 2024)
"c_elegans_uniprot_processing"
#' MEROPS known protease cleavage sites and P4–P4' windows
#'
#' Tabular cleavage events joined to UniProt substrate sequences with P4…P4' windows.
#'
#' @format A data frame with one row per cleavage site.
#' \describe{
#'   \item{protease_id}{MEROPS protease identifier}
#'   \item{uniprot_acc}{UniProt accession of the substrate}
#'   \item{position}{1-based index of the P1' residue in the substrate sequence (MEROPS convention)}
#'   \item{window8}{Eight-residue window P4…P4' around the cleavage. Near termini the window is padded with X}
#'   \item{evidence_type}{MEROPS evidence type; may be NA}
#'   \item{physio_context}{Physiological context; may be NA}
#'   \item{method}{Experimental method; may be NA}
#'   \item{reference}{Primary reference; may be NA}
#' }
#' @details
#' Windows are constructed from the substrate sequence with non-standard letters mapped to X.
#'
#' @source MEROPS database, current_release/database_files (cleavage.txt, substrate.txt);
#'   ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/database_files/ (accessed September 2025).
"merops_sites"
#' Position-specific scoring matrices per MEROPS protease
#'
#' PSSMs built from merops_sites window8 strings, using global amino-acid background
#' and a BLOSUM62-smoothed pseudocount prior. Scores are log2-odds.
#'
#' @format A data frame with:
#' \describe{
#'   \item{protease_id}{MEROPS protease identifier}
#'   \item{n_sites}{Number of cleavage sites used for the model}
#'   \item{pssm}{20×8 numeric matrix; rows are amino acids (A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V),
#'               columns are positions P4,P3,P2,P1,P1',P2',P3',P4'}
#'   \item{bg}{Length-20 numeric vector of global amino-acid background frequencies (same row order as pssm)}
#' }
#' @details
#' PSSM probabilities are computed per column with pseudocount pc=1 mixed from a softened
#' BLOSUM62 prior (temperature=3), then converted to log2 odds vs the global background.
#'
#' @source Derived from \code{merops_sites}. Built by \code{scr/merops_build.R}.
"merops_pssm"
#' MEROPS protease ID to name mapping
#'
#' Two-column lookup mapping MEROPS internal protease identifiers to human-readable protease names.
#'
#' @format A data frame with:
#' \describe{
#'   \item{merops_protease_id}{MEROPS protease identifier (e.g., C14.003)}
#'   \item{protease_name}{Protease name}
#' }
#'
#' @details
#' Derived from the MEROPS uniprot.txt mapping file and deduplicated.
#'
#' @source MEROPS database, current_release/database_files/uniprot.txt;
#'   ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/database_files/ (accessed September 2025).
"merops_accession_to_protease"
