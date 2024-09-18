#' Annotate peptides and associated cleavage events
#'
#' @title annotate_neo_termini
#' @param peptides_df peptide data frame as obtained from one of the adapter functions. Must contain columns nterm_modif_peptide, nterm_modif, peptide, protein and one quantitative column per sample.
#' @param fasta_location location of the fasta file containing the protein sequences identified in the sample, in uniprot format.
#' @param sense direction of the peptide cleavage. One of "N" or "C".
#' @param specificity amino acid specificity of the cleavage site. Example: "K|R" for trypsin.
#' @param organism organism of the fasta file. Default is "mouse". Other options are "human", "mendicato_trucantula", "rhizobium_melitoli", "pig", "human_iso".
#' @param distinct logical. If TRUE, keep only one peptide sequence per feature. Default is TRUE.
#'
#' @format A data frame with at least 28 columns
#' \describe{
#'  \item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
#'  \item{nterm_modif}{N-terminal modification}
#'  \item{peptide}{Peptide sequence}
#'  \item{protein}{Protein ID based on Uniprot Accession}
#'  \item{protein_length}{Protein length}
#'  \item{peptide_start}{Start position of the peptide in the protein sequence}
#'  \item{peptide_end}{End position of the peptide in the protein sequence}
#'  \item{last_aa}{Last amino acid of the peptide}
#'  \item{first_aa}{First amino acid of the peptide}
#'  \item{aa_after}{Amino acid after the peptide}
#'  \item{aa_before}{Amino acid before the peptide}
#'  \item{sense}{Direction of the peptide cleavage ("C" means after C-termini, like trypsin)}
#'  \item{specificity}{Specificity of the cleavage sites of the peptide ("specific", "semi_Nterm", or "semi_Cterm")}
#'  \item{five_res_before}{Five amino acids before the cleavage site}
#'  \item{five_res_after}{Five amino acids after the cleavage site}
#'  \item{cleavage_site}{Five amino acids before and after the cleavage site, separated by "|"}
#'  \item{cleavage_sequence}{Ten amino acids around the cleavage site}
#'  \item{p1_position}{Position of the P1 residue around the cleavage site}
#'  \item{p1_prime_position}{Position of the P1' residue around the cleavage site}
#'  \item{p1_residue}{P1 residue (amino acid)}
#'  \item{p1_prime_residue}{P1' residue (amino acid)}
#'  \item{p1_position_percentage}{Position of the P1 residue in proportion of the protein length}
#'  \item{met_clipping}{Logical. TRUE if the peptide is a methionine clipping}
#'  \item{matches_p1_prime}{Logical. TRUE if the peptide matches the P1' position of a processing feature}
#'  \item{uniprot_processing_type}{Type of processing feature from Uniprot API. One of "INIT_MET", "PROPEP", 
#'                                 "SIGNAL", "TRANSIT", "CHAIN". 
#'                                 "non_canonical" indicates that the potential cleavage location 
#'                                 has not been annotated in UniProt processing features.
#'                                 "not_canonical_no_procc_annot" indicates that the potential cleavage site
#'                                 maps to a Uniprot protein ID that does not have processing features information.}
#'  \item{processing_annotation_start}{Start position of the processing feature in the protein sequence}
#'  \item{processing_annotation_end}{End position of the processing feature in the protein sequence}
#'  \item{protein_sequence}{Protein sequence}
#'  \item{sample_columns}{Scaled abundances of annotated features per sample}
#' }
#'
#' @importFrom dplyr select distinct mutate filter left_join arrange group_by summarize ungroup relocate
#' @importFrom stringr str_extract str_locate str_sub str_length
#' @importFrom seqinr read.fasta
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#'
#' @description
#' Annotate peptides and associated cleavage events.
#'
#'
#' @export
#'
#' @author Miguel Cosenza-Contreras

annotate_neo_termini <- function(
  peptides_df,
  fasta_location,
  sense,
  specificity,
  organism,
  distinct = TRUE){

  require(dplyr)
  require(stringr)
  require(seqinr)
  require(purrr)
  require(tidyr)
  require(magrittr)
  
  data("unimod_id_to_name_mapping", package = "TermineR")

  peptide2protein <- peptides_df %>%
                    dplyr::select(peptide, protein, nterm_modif) %>%
                    dplyr::distinct()

  fasta_2 <- seqinr::read.fasta(
        file = fasta_location,
        seqtype = "AA",
        as.string = TRUE,
        set.attributes = FALSE
    )

  fasta_names <- str_extract(
    names(fasta_2),
    "(?<=\\|)(.*?)(?=\\|)"
    )

  fasta_df <- tibble(
        protein = fasta_names,
        protein_sequence = unlist(fasta_2)
    )
  
  expected_modifications <- c("TMT",
                              "Dimethyl",
                              "Acetyl",
                              str_subset(unimod_id_to_name_mapping$name,
                                         "Dimethyl"),
                              str_subset(unimod_id_to_name_mapping$name,
                                         "Acetyl"),
                              str_subset(unimod_id_to_name_mapping$name,
                                         "TMT"))

prot2pept2fasta <- left_join(
    peptide2protein,
    fasta_df,
    by = c("protein" = "protein")
  ) %>%
  mutate(
    position = str_locate(protein_sequence, peptide),
    start_position = position[, 1],
    end_position = position[, 2],
    last_aa = str_sub(peptide, -1, -1),
    first_aa = str_sub(peptide, 1, 1),
    aa_after = str_sub(protein_sequence, end_position + 1, end_position + 1),
    aa_before = str_sub(protein_sequence, start_position - 1, start_position - 1),
    sense = sense,
    protein_length = str_length(protein_sequence)
  ) %>%
  mutate(
    # substitute "" empty strings with NA for aa_after and aa_before
    aa_after = if_else(aa_after == "", NA_character_, aa_after),
    aa_before = if_else(aa_before == "", NA_character_, aa_before)
  ) %>%
  dplyr::select(-all_of(c("position"))) %>%
  # annotate specificity
  mutate(
    specificity = case_when(
      sense == "C" & str_detect(last_aa, specificity) & str_detect(aa_before, specificity, negate = TRUE) & !is.na(aa_before) ~ "semi_Nterm",
      sense == "C" & str_detect(last_aa, specificity, negate = TRUE) & str_detect(aa_before, specificity) & !is.na(aa_after) ~ "semi_Cterm",
      sense == "N" & str_detect(first_aa, specificity) & str_detect(aa_after, specificity, negate = TRUE) & !is.na(aa_after) ~ "semi_Cterm",
      sense == "N" & str_detect(first_aa, specificity, negate = TRUE) & str_detect(aa_after, specificity) & !is.na(aa_before) ~ "semi_Nterm",
      sense == "C" & str_detect(last_aa, specificity, negate = TRUE) & str_detect(aa_before, specificity, negate = TRUE) & !is.na(aa_before) & !is.na(aa_after) & aa_before == "M" ~ "semi_Cterm",
      sense == "C" & str_detect(aa_before, specificity, negate = TRUE) & !is.na(aa_before) & is.na(aa_after) ~ "semi_Nterm",
      sense == "C" & str_detect(last_aa, specificity, negate = TRUE) & !is.na(aa_after) & is.na(aa_before) ~ "semi_Cterm",
      sense == "C" & str_detect(last_aa, specificity, negate = TRUE) & str_detect(aa_before, specificity, negate = TRUE) & !is.na(aa_before) & !is.na(aa_after) & aa_before != "M" ~ "unspecific",
      sense == "N" & str_detect(first_aa, specificity, negate = TRUE) & str_detect(aa_after, specificity, negate = TRUE) & !is.na(aa_before) & !is.na(aa_after) & aa_before != "M"~ "unspecific",
      TRUE ~ "specific"
  )
  ) %>%
  mutate(
    five_res_before = case_when(
        specificity == "semi_Nterm" & start_position - 5 <= 1 ~ str_sub(protein_sequence, 1, start_position - 1),
        specificity == "semi_Nterm" & start_position - 5 > 1 ~ str_sub(protein_sequence, start_position - 5, start_position - 1),
        #specificity == "semi_Cterm" & end_position + 5 >= str_length(protein_sequence) ~ str_sub(protein_sequence, end_position - 4, str_length(protein_sequence)),
        #specificity == "semi_Cterm" & end_position + 5 < str_length(protein_sequence) ~ str_sub(protein_sequence, end_position - 4, end_position),
        specificity == "semi_Cterm" ~ str_sub(protein_sequence, end_position - 4, end_position),
        specificity == "specific" & nterm_modif %in% expected_modifications & start_position - 5 <= 1 ~ str_sub(protein_sequence, 1, start_position - 1),
        specificity == "specific" & nterm_modif %in% expected_modifications & start_position - 5 > 1 ~ str_sub(protein_sequence, start_position - 5, start_position - 1),
        specificity == "specific" & !nterm_modif %in% expected_modifications & start_position - 5 <= 1 ~ str_sub(protein_sequence, 1, start_position - 1),
        specificity == "specific" & !nterm_modif %in% expected_modifications & start_position - 5 > 1 ~ str_sub(protein_sequence, start_position - 5, start_position - 1)
    ),
    five_res_after = case_when(
        specificity == "semi_Nterm" & start_position + 4 < str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, start_position + 4),
        specificity == "semi_Nterm" & start_position + 4 >= str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, str_length(protein_sequence)),
        specificity == "semi_Cterm" & end_position + 4 < str_length(protein_sequence) ~ str_sub(protein_sequence, end_position + 1, end_position + 5),
        specificity == "semi_Cterm" & end_position + 4 >= str_length(protein_sequence) ~ str_sub(protein_sequence, end_position + 1, str_length(protein_sequence)),
        specificity == "specific" & nterm_modif %in% expected_modifications & start_position + 4 < str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, start_position + 4),
        specificity == "specific" & nterm_modif %in% expected_modifications & start_position + 4 >= str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, str_length(protein_sequence)),
        specificity == "specific" & !nterm_modif %in% expected_modifications & start_position + 4 < str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, start_position + 4),
        specificity == "specific" & !nterm_modif %in% expected_modifications & start_position + 4 >= str_length(protein_sequence) ~ str_sub(protein_sequence, start_position, str_length(protein_sequence))
    )
    ) %>%
  mutate(
    five_res_after = str_pad(five_res_after, 5, side = "right", pad = "X"),
    five_res_before = str_pad(five_res_before, 5, side = "left", pad = "X")
  ) %>%
  mutate(
    cleavage_site = paste0(
        five_res_before,
        " | ",
        five_res_after
      ),
    cleavage_sequence = paste0(
        five_res_before,
        five_res_after
      ),
    p1_position = case_when(
      specificity == "semi_Nterm" ~ start_position - 1,
      specificity == "semi_Cterm" ~ end_position,
      specificity == "specific" ~ start_position - 1,
      TRUE ~ NA
    ),
    p1_prime_position = case_when(
      specificity == "semi_Nterm" ~ start_position,
      specificity == "semi_Cterm" ~ end_position + 1,
      specificity == "specific" ~ start_position,
      TRUE ~ NA
    )
  ) %>%
  mutate(
    p1_residue = case_when(
      specificity != "specific" ~ str_sub(five_res_before, 5, 5),
      specificity == "specific" & sense == "C" ~ aa_before,
      specificity == "specific" & sense == "N" ~ last_aa,
      p1_position == 0 ~ "N-term"
      TRUE ~ NA),
    p1_prime_residue = case_when(
      specificity != "specific" ~ str_sub(five_res_after, 1, 1),
      specificity == "specific" & sense == "C" ~ first_aa,
      specificity == "specific" & sense == "N" ~ aa_after,
      p1_prime_position == protein_length ~ "C-term",
      TRUE ~ NA),
    p1_position_percentage = p1_position / str_length(protein_sequence) * 100,
    met_clipping = case_when(
      specificity != "specific" & start_position <= 3 & (str_sub(five_res_before, 5, 5) == "M" | str_sub(five_res_before, 4, 4) == "M") ~ TRUE,
      TRUE ~ FALSE
    )
  )

annotated_df_w_quant <- left_join(
    peptides_df,
    prot2pept2fasta,
    by = c("peptide" = "peptide", "protein" = "protein", "nterm_modif" = "nterm_modif"),
    multiple = "first"
  )

  # annotate cleavage sites based on uniprot processing database

protein_nterm_modif <- annotated_df_w_quant %>%
  dplyr::select(
    protein,
    peptide,
    nterm_modif,
    specificity,
    last_aa,
    aa_before,
    peptide_start = start_position,
    peptide_end = end_position,
    met_clipping,
    p1_position,
    p1_prime_position,
    p1_residue,
    p1_prime_residue,
    protein_length) %>%
  dplyr::filter(str_detect(protein,
                           pattern = "Biognosys",
                           negate = TRUE),
                nterm_modif != "n")

protein_semi_free <- annotated_df_w_quant %>%
  dplyr::select(
  protein,
  peptide,
  nterm_modif,
  specificity,
  last_aa, aa_before,
  peptide_start = start_position,
  peptide_end = end_position,
  met_clipping,
  p1_position,
  p1_prime_position,
  p1_residue,
  p1_prime_residue,
  protein_length) %>%
  dplyr::filter(str_detect(protein,
                           pattern = "Biognosys",
                           negate = TRUE),
                specificity %in% c("semi_Nterm", "semi_Cterm"),
                nterm_modif %in% c("n"))

protein_nter <- bind_rows(
  protein_nterm_modif,
  protein_semi_free
  )

# vector of interesting feature types as annotated by the uniprot API
mol_processing_feat <- c(
  "INIT_MET",
  "PROPEP",
  "SIGNAL",
  "TRANSIT",
  "CHAIN"
  )

expected_organisms <- c(
  "mouse",
  "human",
  "mendicato_trucantula",
  "rhizobium_melitoli",
  "pig",
  "human_iso"
)

if(organism == "mouse"){

  data("mouse_uniprot_processing", package = "TermineR")

  uniprot_processing <- mouse_uniprot_processing

  } else if (organism == "human"){
  
  data("human_uniprot_processing", package = "TermineR")

  uniprot_processing <- human_uniprot_processing
  
  } else if(organism == "mendicato_trucantula"){

  data("mendicato_trucantula_uniprot_processing", package = "TermineR")

  uniprot_processing <- mendicato_trucantula_uniprot_processing

  } else if(organism == "rhizobium_melitoli"){

  data("rhizobium_melitoli_uniprot_processing", package = "TermineR")

  uniprot_processing <- rhizobium_melitoli_uniprot_processing

  } else if(organism == "pig"){

  data("pig_uniprot_processing", package = "TermineR")

  uniprot_processing <- pig_uniprot_processing

  } else if(organism == "human_iso"){
    
  data("human_and_isoforms_uniprot_processing", package = "TermineR")
    
  uniprot_processing <- human_and_isoforms_uniprot_processing
  
  } else if(organism %in% expected_organisms == FALSE){
    
  stop("Organism not found in our list of annotations. Check documentation for list of supported organisms.")
  
  }

df_mol_proc_feat <- uniprot_processing %>%
  dplyr::filter(
    type %in% mol_processing_feat) %>% #, # keep only interesting features
                #!is.na(length)) %>% # exclude features with missing values
  dplyr::rename(protein = accession)  # change column name

nter_pepts_n_feat <- left_join(
  protein_nter,
  df_mol_proc_feat,
  by = "protein",
  relationship = "many-to-many")

# match semi-specific cleavage position of peptide vs end position of annotated processing site

categ_canon_annot <- nter_pepts_n_feat %>%
    mutate(
      matches_p1_prime = case_when(
        abs(as.numeric(p1_prime_position) - end) < 4 ~ TRUE,
        TRUE ~ FALSE
        )
        ) %>%
  #dplyr::filter(!is.na(type)) %>% # eliminate proteins with no processing features annotated
  mutate(
    processing_type = case_when(
      matches_p1_prime == TRUE & type == "CHAIN" ~ "CHAIN",
      matches_p1_prime == TRUE & type == "INIT_MET" ~ "INIT_MET",
      matches_p1_prime == TRUE & type == "SIGNAL" ~ "SIGNAL",
      matches_p1_prime == TRUE & type == "TRANSIT" ~ "TRANSIT",
      matches_p1_prime == TRUE & type == "PROPEP" ~ "PROPEP",
      matches_p1_prime == TRUE & type == "PEPTIDE" ~ "PEPTIDE",
      matches_p1_prime == FALSE ~ "not_canonical",
      is.na(matches_p1_prime) ~ "not_canonical_no_procc_annot",
      TRUE ~ "not_canonical_no_procc_annot"
    )
  ) %>%
  # select interesting columns
  dplyr::select(
    peptide,
    protein,
    matches_p1_prime,
    processing_type,
    p1_position,
    p1_prime_position,
    p1_prime_residue,
    p1_residue,
    peptide_start,
    peptide_end,
    processed_fragment_start = begin,
    processed_fragment_end = end,
    nterm_modif,
    specificity,
    met_clipping) %>%
  dplyr::filter(
    processing_type != "CHAIN"
  )

# filter to keep one peptide sequence per feature
pept_wmatch <- categ_canon_annot %>%
    dplyr::filter(
      !processing_type %in% c("not_canonical", "not_canonical_no_procc_annot")
    )

  if(distinct == TRUE){

  pept_wmatch <- pept_wmatch %>%
      distinct(
        peptide,
        protein,
        matches_p1_prime,
        .keep_all = TRUE
      )

    } else {

    pept_wmatch <- pept_wmatch

    }

pept_womatch <- categ_canon_annot %>%
    dplyr::filter(!peptide %in% pept_wmatch$peptide) %>%
    distinct(peptide, protein,
             matches_p1_prime, processing_type,
             .keep_all = TRUE)

categ2_pept_canannot <- bind_rows(pept_wmatch,
                                  pept_womatch) %>%
    mutate(is_duplicated = duplicated(peptide)) %>%
    # adjust the definition of INIT_MET based on M at P1
    mutate(
      processing_type = case_when(
        processing_type == "not_canonical" & p1_residue == "M" & peptide_start == 2 ~ "INIT_MET_not_canonical",
        processing_type == "not_canonical" & p1_prime_position == 1 ~ "Intact_ORF",
        processing_type == "not_canonical_no_procc_annot" & p1_residue == "M" & peptide_start == 2 ~ "INIT_MET_not_canonical",
        processing_type == "not_canonical_no_procc_annot" & p1_prime_position == 1 ~ "Intact_ORF",
        TRUE ~ processing_type
      ),
    ) %>%
    dplyr::select(
      peptide,
      protein,
      matches_p1_prime,
      uniprot_processing_type = processing_type,
      processing_annotation_start = processed_fragment_start,
      processing_annotation_end = processed_fragment_end
    )

final_annotated_df <- left_join(
    annotated_df_w_quant %>% dplyr::rename(
      peptide_start = start_position, 
      peptide_end = end_position),
    categ2_pept_canannot,
    relationship = "many-to-many"
    ) %>%
  mutate(
    uniprot_processing_type = case_when(
      is.na(matches_p1_prime) ~ "not_canonical_no_procc_annot",
      TRUE ~ uniprot_processing_type
    )
  ) %>%
  dplyr::relocate(
    nterm_modif_peptide,
    nterm_modif,
    peptide,
    protein,
    protein_length,
    peptide_start,
    peptide_end,
    last_aa,
    first_aa,
    aa_after,
    aa_before,
    sense,
    specificity,
    five_res_before,
    five_res_after,
    cleavage_site,
    cleavage_sequence,
    p1_position,
    p1_prime_position,
    p1_residue,
    p1_prime_residue,
    p1_position_percentage,
    met_clipping,
    matches_p1_prime,
    uniprot_processing_type,
    processing_annotation_start,
    processing_annotation_end,
    protein_sequence
  )

  return(final_annotated_df)

}

