#' Normalize peptide abundances against their respective protein abundances.
#'
#' @title peptide2protein_normalization
#' @param peptides A data frame with the raw (not log2 transformed) peptide abundances per sample. It should contain the columns 'nterm_modif_peptide', 'protein' and the sample names as columns.
#' @param annot A data frame with the annotation of the samples. It should contain at least a 'sample' column, matching the column names for the samples in the peptide data frame.
#' @param peptide_annot A data frame with the annotation of the peptides. It should contain at least a 'nterm_modif_peptide' column, matching the 'nterm_modif_peptide' column in the peptide data frame, and the 'specificity' column.
#' @param summarize_by_specificity A logical value indicating if the protein abundances should be summarized based on specific peptides. Default is TRUE.
#'
#' @description This experimental function uses summed fully specific peptide abundances (fully tryptic for a trypsin search) as a proxy for parent-protein abundance when summarize_by_specificity = TRUE. Each peptide abundance is divided by its corresponding protein proxy, then log2 transformed and median centered. The output is a protein-adjusted peptide-abundance proxy, not an absolute molecular fraction or a direct measure of protease activity. No external standards, reference samples, control groups, or reference peptide/protein values are required.
#'
#' @importFrom dplyr select filter pull mutate group_by summarise_if rename_at arrange slice relocate
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom stringr str_remove
#'
#' @return A list with the following elements
#' \describe{
#'  \item{protein_normalized_pepts_scaled}{Protein-adjusted peptide-abundance proxies, log2 transformed and median centered}
#'  \item{protein_normalized_pepts_abundance}{Peptide/protein abundance ratios representing protein-adjusted peptide-abundance proxies, before log2 transformation}
#'  \item{summarized_protein_abundance}{Summarized protein abundances based on peptide matrix}
#'  \item{summarized_protein_abundance_scaled}{Summarized protein abundances based on peptide matrix, scaled}
#'  \item{summarize_by_specificity}{Logical flag indicating if protein abundances were summarized by fully specific peptides}
#' }
#'
#' @export
#' @author Miguel Cosenza-Contreras
peptide2protein_normalization <- function(peptides,
                                          annot,
                                          peptide_annot,
                                          summarize_by_specificity = TRUE) {

  require(dplyr)
  require(purrr)
  require(tidyr)
  require(tibble)
  require(magrittr)

  # pre-process annotation file --------------------------

  # we need to modify the sample names in the annotation file to match the
  # names of the of the TMT intensities after loading into R.

  psm_cols_trim <- colnames(peptides) %>% # get column names of psm file
    str_remove("\\..*") # eliminate sufix

  psm_cols <- colnames(peptides) # get column names of psm file wo eliminate sufix

  # create df mapping trimmed sample names and the ones with sufix
  trimed2original_sample <- tibble(sample = psm_cols_trim,
                                   sample_trim = psm_cols)

  # correct annotation file to match column names for TMT channel intensities
  # in the psm file
  annot_2 <- annot %>%
    left_join(., trimed2original_sample) %>% # merge annotation file w untrimmed
    dplyr::slice(1:nrow(annot)) %>% # keep only sample names present in columns of psm file
    # elimitate 'empty' samples
    filter(sample != "empty")

  # generate protein-quant data frame from peptide abundances ----

  ## first, we want to be able to differentiate between specific or semi-specific peptides
  ## ideally, we would want to normalize peptide abundances against protein
  ## abundances summarized from fully specific peptides.
  ## This is because we assume that protein abundance based on specific peptides
  ## better resembles overall abundance of the protein, without the effect of
  ## proteolytic activity.

  if (summarize_by_specificity == TRUE){

    specific_peptides <- peptide_annot %>%
      filter(specificity == "specific") %>%
      pull(nterm_modif_peptide)

    prots_q <- peptides %>%
      # keep only unique peptides/PSMs
      filter(nterm_modif_peptide %in% specific_peptides) %>%
      dplyr::select(protein,
                    nterm_modif_peptide,
                    annot_2$sample_trim) %>%
      group_by(protein) %>%
      summarise_if(is.numeric,
                   sum,
                   na.rm = TRUE) %>%
      rename_at(.vars = vars(annot_2$sample_trim),
                .funs = function(x) paste0(x, "_prot"))

  } else if (summarize_by_specificity == FALSE){

    prots_q <- peptides %>%
      dplyr::select(protein,
                    nterm_modif_peptide,
                    annot_2$sample_trim) %>%
      group_by(protein) %>%
      summarise_if(is.numeric,
                   sum,
                   na.rm = TRUE) %>%
      rename_at(.vars = vars(annot_2$sample_trim),
                .funs = function(x) paste0(x, "_prot"))

  }

  # processing files -----

  ## select quant columns from peptide data frame and add suffix

  peptides_q <- peptides %>%
    dplyr::select(nterm_modif_peptide, protein,
                  annot_2$sample_trim) %>%
    rename_at(.vars = vars(annot_2$sample_trim),
              .funs = function(x) paste0(x, "_peptide")) # add prot suffix to columns

  ## merge peptide file with protein file

  psm_n_prots <- left_join(peptides_q, prots_q) %>%
    dplyr::relocate(nterm_modif_peptide)

  # function to get the ratio of intensity of Peptide/Protein  -----

  ## this function calculates peptide abundance relative to a parent-protein proxy
  ## calculated from fully specific peptides when summarize_by_specificity is TRUE

  pept2prot_ratios <- function(col){

    df_q <- psm_n_prots %>%
      # select id column and columns that match sample names (for map function below)
      dplyr::select(nterm_modif_peptide, starts_with(col)) %>%
      # reorganize the column of intensities, so first we have peptide abundances and then proteins
      dplyr::relocate(nterm_modif_peptide, ends_with("peptide"), ends_with("prot"))

    df_q_names <- colnames(df_q) # extract column names

    df_q_rat <- df_q %>%
      # create a peptide/protein ratio for each sample; invalid protein denominators yield NA
      mutate({{col}} := ifelse(
        is.finite(.data[[df_q_names[3]]]) & .data[[df_q_names[3]]] > 0,
        .data[[df_q_names[2]]] / .data[[df_q_names[3]]],
        NA_real_
      )) %>%
      # retain the ratio as a protein-adjusted peptide-abundance proxy
      mutate("fraction_int_peptide2prot_{col}" := .data[[{{col}}]])
    df_q_rat <- df_q_rat %>%
      # keep only the columns with peptide/protein abundance ratios
      dplyr::select(-c(nterm_modif_peptide, ends_with("peptide"), ends_with("prot")))
    return(df_q_rat)

  }

  ## PSM intensity/protein ration calculation per sample ----

  list_pept2prot_rat <- purrr::map(.x = annot_2$sample,
                                   .f = pept2prot_ratios)

  pept2prot_norm1 <- bind_cols(list_pept2prot_rat) %>%
    mutate(nterm_modif_peptide = psm_n_prots$nterm_modif_peptide) %>%
    dplyr::relocate(nterm_modif_peptide)

  pept2prot_norm_ratio <- pept2prot_norm1 %>%
    dplyr::select(nterm_modif_peptide, starts_with("fraction_int_"))

  ## normalizations ---------------------------------------------------

  ### log2 and median centering of peptide/protein abundance ratios -----

  mat_ratios2 <- pept2prot_norm_ratio %>%
    column_to_rownames("nterm_modif_peptide") %>%
    as.matrix()

  mat_ratios2[mat_ratios2 == 0] <- NA

  log2_mat_rat2 <- mutate_all(as.data.frame(mat_ratios2),
                              log2)

  scaled_mat_rat2 <- scale(log2_mat_rat2,
                           scale = F,
                           center = apply(log2_mat_rat2, 2, median,
                                          na.rm = TRUE) - median(as.matrix(log2_mat_rat2),
                                                                 na.rm = TRUE)) %>%
    # log2 values may be negative; retain their signs
    na.omit(.) %>% # remove rows with NA values
    as.data.frame(.) %>%
    rownames_to_column("nterm_modif_peptide")

  ### log2 and median centering of summarized protein abundances -----

  mat_ratios3 <- prots_q %>%
    column_to_rownames("protein") %>%
    as.matrix()

  mat_ratios3[!is.finite(mat_ratios3) | mat_ratios3 <= 0] <- NA

  log2_mat_rat3 <- mutate_all(as.data.frame(mat_ratios3),
                              log2)

  scaled_mat_rat3 <- scale(log2_mat_rat3,
                           scale = F,
                           center = apply(log2_mat_rat3,
                                          2,
                                          median,
                                          na.rm = TRUE) - median(as.matrix(log2_mat_rat3),
                                                                 na.rm = TRUE)) %>%
    # log2 values may be negative; retain their signs
    as.data.frame() %>%
    rownames_to_column("protein")


  # create list of results ----

  protein_normalized_pepts <- list(
    # protein-adjusted peptide-abundance proxies, log2 transformed and median centered
    protein_normalized_pepts_scaled = scaled_mat_rat2,
    # protein-adjusted peptide-abundance proxies as peptide/protein ratios before log2 transformation
    protein_normalized_pepts_abundance = pept2prot_norm_ratio,
    # summarized protein abundances based on peptide matrix
    summarized_protein_abundance = prots_q,
    # summarized protein abundances based on peptide matrix, scaled
    summarized_protein_abundance_scaled = scaled_mat_rat3,
    # object showing if the protein abundances were summarized by tryptic peptides
    summarize_by_specificity = summarize_by_specificity)

  return(protein_normalized_pepts)
}

