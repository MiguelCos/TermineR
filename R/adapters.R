#' Process psm.tsv files from FragPipe search results for further analysis with TermineR
#'
#' @title fragpipe_adapter
#' @param parent_dir directory of FragPipe search results
#' @param ref_sample define the reference channel in multi-mixture TMT experiments. Default is NULL (1 mixture). Virtual reference is not available.
#' @param grouping_var define the grouping variable for MAD scaling. Default is "nterm_modif_peptide" for normalization at the feature level of "Peptide + N-terminal modification".
#' @param min_purity minimum purity of the PSMs to be considered. Default is 0.5.
#' @param tmt_delta define the TMT reporter ion mass delta. Default is "229" for TMT10/11plex. "304" for TMT16plex.
#'
#' @format A data frame with at least 4 columns
#' \describe{
#'  \item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
#'  \item{nterm_modif}{N-terminal modification}
#'  \item{peptide}{Peptide sequence}
#'  \item{protein}{Protein ID based on Uniprot Accession}
#' }
#'
#'
#' @description
#' Process psm.tsv files from FragPipe search results for further analysis with TermineR
#'
#' @importFrom dplyr mutate filter select bind_rows left_join group_by summarize ungroup arrange distinct slice_max
#' @importFrom tidyr pivot_wider separate
#' @importFrom stringr str_detect str_remove
#' @importFrom readr read_tsv
#' @importFrom purrr map map2
#' @importFrom here here
#' @importFrom magrittr %>%
#'
#' @export
#' @author Miguel Cosenza-Contreras
fragpipe_adapter <- function(parent_dir,
                             ref_sample = "norm",
                             grouping_var = "nterm_modif_peptide",
                             min_purity = 0.5,
                             tmt_delta = "229") {

  # required packages

  require(here)
  require(readr)
  require(purrr)
  require(dplyr)
  require(stringr)
  require(tidyr)

  # required functions ------

  ## function to read the annotation files

  read_annotation_files <- function(x){

      require(readr)
      require(dplyr)

      # read the annotation file
      annot <- read_delim(file = x,
                          delim = " ",
                          col_names = FALSE) %>%
      dplyr::rename(channel =  X1,
                    sample_name = X2)

      return(annot)

  }

  # function to filter the psm files from fragpipe
  # based on purity and sum of reporter ion intensities

  filter_psms <- function(x,
                          min_purity = 0.5,
                          tmt_delta = "229",
                          ref_sample = "norm"){

    require(dplyr)
    require(stringr)

    # tmt delta: 229 for TMT10/11plex

    # check if there is an 'empty' column in the psm dataframe and the annotation

    if(any(str_detect(colnames(x$psm), "empty"))){

      x$psm <- x$psm %>%
        dplyr::select(-matches("empty"))

    }

    if("empty" %in% x$annot$sample_name){

      x$annot <- x$annot %>%
        dplyr::filter(sample_name != "empty")

    }

    x$psm <- x$psm %>%
      # in all_of sample columns, substitute zero values with NA
      mutate(across(all_of(x$annot$sample_name),
                    ~case_when(. == 0 ~ NA,
                               TRUE ~ .)))

    # keep interesting columns
    sample_cols <- x$annot$sample_name

    interest_cols <- c("Peptide",
                       "Modified Peptide",
                       "PeptideProphet Probability",
                       "Intensity",
                       "Assigned Modifications",
                       "Purity",
                       "Protein ID",
                       "Gene")

    cols_to_keep <- c(sample_cols, interest_cols)

    # select the psm dataframe to keep only the selected columns
    filtered_psm <- x$psm %>%
      dplyr::select(all_of(cols_to_keep)) %>%
      mutate(sum_reporter_ions = rowSums(dplyr::select(., all_of(sample_cols)), na.rm = TRUE)) %>%
      # create new peptide identification based on N-terminal modifications + peptide sequence
      mutate(
      # annotate N-termini
      nterminal_modification = case_when(
        str_detect(`Assigned Modifications`, "N-term\\(304.207[0-9]\\)") ~ "TMT",
        str_detect(`Assigned Modifications`, "N-term\\(229.162[0-9]\\)") ~ "TMT",
        str_detect(`Assigned Modifications`, "N-term\\(42.010[0-9]\\)") ~ "Acetyl",
        str_detect(`Assigned Modifications`, "N-term\\(28.031[0-9]\\)") ~ "Dimethyl",
        str_detect(`Assigned Modifications`, "N-term\\(36.075[0-9]\\)") ~ "Dimethyl",
        TRUE ~ "n"
        )
      ) %>%
    mutate(
      nterm_modif_peptide = paste(
        nterminal_modification,
        Peptide,
        sep = "_"
      )
    )

    # filter by criteria
    filtered_psm <- filtered_psm %>%
      #filter(rowSums(dplyr::select(., all_of(sample_cols)) == 0) == 0) %>%
      filter(str_detect(string = `Assigned Modifications`,
                        pattern = tmt_delta),
             Purity > min_purity,
             sum_reporter_ions > 0) %>%
    # if several PSMs to the same nterminal modified peptide, select one with best Purity
    group_by(nterm_modif_peptide) %>%
    slice_max(n = 1,
              order_by = c(Purity)) %>%
    # if a modified peptide was identified with several PSMs and same high purity,
    # then select the PSMs with highest summed TMT intensity
    slice_max(n = 1,
              order_by = sum_reporter_ions) %>%
    ungroup()

    if(!is.null(ref_sample)){

    filtered_psm <- filtered_psm %>%
      filter(.data[[ref_sample]] > 0)
    }

    annot <- x$annot

    # Return the list with the filtered psm dataframe and the original annot dataframe
    return(list(psm = filtered_psm,
                annot = annot))

  }

  # mad scaling function
  # this function will correct for reference channel within each mixture
  # and summarize abundances by modified peptide, peptide, protein, and gene

  mad_scaling <- function(x,
                          ref_sample = "sc",
                          grouping_var = "Protein ID",
                          mixture_col) {

    require(dplyr)
    require(tidyr)

    # define interesting columns
    interest_cols <- c("Peptide",
                       "Modified Peptide",
                       "PeptideProphet Probability",
                       "Intensity",
                       "Assigned Modifications",
                       "Purity",
                       "Protein ID",
                       "Gene")

    # keep interesting columns
    sample_cols <- x$annot$sample_name

    cols_to_keep <- c(sample_cols, interest_cols)

    if(!is.null(ref_sample)){

    # select the psm dataframe to keep only the selected columns
    scaled_ratios <- x$psm %>%
    mutate(across(all_of(sample_cols),
                  ~log2(. + 1) - log2(.data[[ref_sample]] + 1),
                  .names = "ratio_{.col}"))
    } else {

    # select the psm dataframe to keep only the selected columns
    scaled_ratios <- x$psm %>%
    # in all_of sample columns, substitute zero values with NA
    mutate(across(all_of(sample_cols),
                    ~case_when(. == 0 ~ NA,
                               TRUE ~ .))) %>%
    mutate(across(all_of(sample_cols),
                  ~log2(. + 1),
                  .names = "ratio_{.col}"))

    }

    # exclude ref_sample column

    if(!is.null(ref_sample)){

    scaled_ratios <- scaled_ratios %>%
      dplyr::select(-starts_with(paste0("ratio_", ref_sample)))

    }

    scaled_ratios <- scaled_ratios %>%
    dplyr::select(all_of(c(grouping_var)), starts_with("ratio_")) %>%

    # transform to long format to facilitate calculations at different levels
    pivot_longer(cols = starts_with("ratio_"),
                 names_to = "sample",
                 values_to = "log2_ratio_2_ref") %>%

    # calculate median ratios per desired level, per sample
    group_by(.data[[grouping_var]], sample) %>%
    summarize(log2_rat2ref_group = median(log2_ratio_2_ref,
                                          na.rm = TRUE),
              .groups = "drop") %>%

    # calculate median ratios per sample
    group_by(sample) %>%
    mutate(Mi = median(log2_rat2ref_group,
                       na.rm = TRUE)) %>%
    ungroup() %>%

    # global median across all samples
    mutate(M0 = median(Mi,
                       na.rm = TRUE)) %>%

    # median centered ratios based on median ratios per sample
    mutate(RCij = log2_rat2ref_group - Mi) %>%

    # calculate median absolute deviation per sample
    group_by(sample) %>%
    mutate(MADi = median(abs(RCij),
                         na.rm = TRUE)) %>%
    ungroup() %>%

    # calculate global median absolute deviation
    mutate(MAD0 = median(MADi,
                         na.rm = TRUE)) %>%
    ungroup() %>%

    # calculate scaled ratios
    mutate(RNij = (RCij / MADi) * MAD0 + M0) %>%
    mutate(mixture = mixture_col)

    return(scaled_ratios)

  }

  # function to calculate the reference intensity across mixtures

  calculate_ref_intensity <- function(x,
                                      grouping_var = "Protein ID",
                                      ref_sample = "sc",
                                      mixture_col) {

    require(dplyr)
    require(tidyr)

    if(!is.null(ref_sample)){

    # keep interesting columns
    sample_cols <- x$annot$sample_name

    # define interesting columns
    interest_cols <- c("Peptide",
                       "Modified Peptide",
                       "PeptideProphet Probability",
                       "Intensity",
                       "Assigned Modifications",
                       "Purity",
                       "Protein ID",
                       "Gene")

    cols_to_keep <- c(sample_cols, interest_cols)

    # select the psm dataframe to keep only the selected columns
    reference_intensity <- x$psm %>%
    dplyr::select(all_of(c(grouping_var, sample_cols, "Intensity"))) %>%
    # select the top 3 most intense PSMs per grouping variable
    slice_max(order_by = Intensity,
              by = .data[[grouping_var]],
              n = 3) %>%

    # calculate the sum of reporter ions
    mutate(sum_reporter_ions = rowSums(dplyr::select(.,
                                                     all_of(sample_cols)),
                                       na.rm = TRUE)) %>%

    # calculate the weighted intensity
    mutate(weighted_intensity = Intensity * (.data[[ref_sample]] / sum_reporter_ions)) %>%

    # group by the grouping variable and calculate the reference intensity per grouping variable in each mixture
    group_by(.data[[grouping_var]]) %>%
    summarise(REFik = sum(weighted_intensity),
              .groups = 'drop') %>%
    # add the mixture column
    mutate(mixture = mixture_col)
    } else {

      reference_intensity <- NULL

    }


  return(reference_intensity)

  }

  # Function to calculate the final abundance values including the refenrece intensity

  calculate_final_abundance <- function(scaled_ratios = scaled_ratios,
                                        reference_intensity = reference_intensity,
                                        ref_sample = NULL,
                                        grouping_var = "nterm_modif_peptide") {

    require(dplyr)
    require(tidyr)

    if(!is.null(ref_sample)){

    # merge the reference intensity dataframes from reference_intensity list
    ref_int_df <- bind_rows(reference_intensity) %>%

    # calculate the overall reference intensity per grouping variable
    group_by(.data[[grouping_var]]) %>%
    mutate(REFi = mean(REFik)) %>%
    ungroup()

    # merge the scaled ratios dataframes for each mixture with the reference intensity dataframe
    final_abundance <- left_join(
      bind_rows(scaled_ratios), # scaled_ratios dataframes
      ref_int_df,
      by = c(grouping_var, "mixture")
    ) %>%
    # calculate the final abundance values
    mutate(ref_normalized_abundance = RNij + log2(REFi + 1)) %>%
     # substitute -Inf or Inf values with NA for RNij
    mutate(RNij = case_when(
      RNij == -Inf ~ NA,
      RNij == Inf ~ NA,
      TRUE ~ RNij
    ))

    } else {

      # merge the scaled ratios dataframes for each mixture with the reference intensity dataframe
    final_abundance <- bind_rows(scaled_ratios) %>% # scaled_ratios dataframes

    # calculate the final abundance values
    mutate(ref_normalized_abundance = RNij) %>%
     # substitute -Inf or Inf values with NA for RNij
    mutate(RNij = case_when(
      RNij == -Inf ~ NA,
      RNij == Inf ~ NA,
      TRUE ~ RNij
    ))

    }

    return(final_abundance)

  }

  # start execution ------

  # load the data

  # Construct the path to the psm.tsv files

  if(!any(str_detect(list.files(parent_dir), "psm.tsv"))){

    # get the folder names of main project diretory
    folders_dir <- list.dirs(parent_dir,
                             full.names = TRUE,
                             recursive = FALSE)

    intern_foldrs <- list.dirs(parent_dir,
                               full.names = FALSE,
                               recursive = FALSE)

    psm_file_path <- paste0(folders_dir, "/psm.tsv")

  } else {

    intern_foldrs <- "mix_1"

    folders_dir <- parent_dir

    psm_file_path <- paste0(parent_dir, "/psm.tsv")

  }

  # get a list of psm.tsv files
  message("Loading psm.tsv files...")
  list_psms <- purrr::map(.x = psm_file_path,
                          .f = read_tsv,
                          .progress = TRUE) %>%
                          suppressMessages()

  # define names of list of psm files
  names(list_psms) <- intern_foldrs

  #### Annotation files

  # define location
  message("Loading annotation files...")
  annot_file_path2 <- paste0(folders_dir,
                            "/annotation.txt")

  # load annotation files
  list_annot <- purrr::map(.x = annot_file_path2,
                           .f = read_annotation_files,
                          .progress = TRUE) %>%
                          suppressMessages()

  names(list_annot) <- intern_foldrs

  #### combine the psm table with corresponding annotation
  combined_list <- map2(list_psms,
                        list_annot,
                        ~list(psm = .x,
                              annot = .y),
                              .progress = TRUE) %>%
                          suppressMessages()

  # filter the PSM files and generate summarization based on
  # Nterminal modification + peptide sequence

  message("Filtering PSM files...")

  psms_cols <- purrr::map(.x = combined_list,
                          .f = filter_psms,
                          ref_sample = ref_sample,
                          tmt_delta = tmt_delta,
                          min_purity = min_purity,
                          .progress = TRUE)

  # calculate the MAD scaling factors

  message("Calculating reference scaling factors and normalizing...")
  scaled_ratios <- purrr::map2(psms_cols,
                               names(psms_cols),
                               ~mad_scaling(.x,
                                            ref_sample = ref_sample,
                                            grouping_var = grouping_var,
                                            mixture_col = .y),
                               .progress = TRUE)

  reference_intensity <- purrr::map2(psms_cols,
                                     names(psms_cols),
                                     ~calculate_ref_intensity(.x,
                                                              ref_sample = ref_sample,
                                                              grouping_var = grouping_var,
                                                              mixture_col = .y))

  message("Calculating final scaled abundances...")
  final_abundance <-  calculate_final_abundance(scaled_ratios = scaled_ratios,
                                                reference_intensity = reference_intensity,
                                                ref_sample = ref_sample,
                                                grouping_var = grouping_var)

  df_final_abundance <- final_abundance %>%
    dplyr::select(all_of(c(grouping_var, "sample", "ref_normalized_abundance"))) %>%
    mutate(
      sample = str_remove(sample, "ratio_"),
    ) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "ref_normalized_abundance"
    )

  peptide2protein_map <- bind_rows(
    list_psms
  ) %>%
  dplyr::select(
    peptide = Peptide,
    protein = `Protein ID`) %>%
  distinct()

  df_mat_abundance <- df_final_abundance %>%
    tidyr::separate(col = nterm_modif_peptide,
                    into = c("nterm_modif", "peptide"),
                    sep = "_",
                    remove = FALSE) %>%
    left_join(peptide2protein_map, by = "peptide") %>%
    dplyr::relocate(protein, .after = peptide) %>%
    # arrange by peptide sequence
    arrange(peptide)

  return(df_mat_abundance)

}
#' Process report file from Spectronaut search results for further analysis with TermineR
#'
#' @title spectronaut_adapter
#'
#' @param path_to_file path to file with Spectronaut search results in tsv format
#' @param proteotypic keep only proteotypic peptides. Default is TRUE.
#'
#' @format A data frame with at least 4 columns
#' \describe{
#'  \item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
#'  \item{nterm_modif}{N-terminal modification. Only Acetyl annotation is available in the current version}
#'  \item{peptide}{Peptide sequence}
#'  \item{protein}{Protein ID based on Uniprot Accession}
#' }
#'
#'
#' @description
#' Process report file from Spectronaut search results for further analysis with TermineR
#' @importFrom dplyr select mutate filter group_by summarise ungroup rename arrange relocate where
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom stringr str_remove str_remove_all str_count str_detect str_trim
#' @importFrom readr read_tsv
#' @importFrom janitor clean_names
#' @importFrom magrittr %>%
#' @importFrom tidyselect where
#'
#' @export
#' @author Miguel Cosenza-Contreras
spectronaut_adapter <- function(
    path_to_file,
    proteotypic = TRUE) {

  # required packages

  require(here)
  require(readr)
  require(dplyr)
  require(stringr)
  require(tidyr)
  require(janitor)

  # required functions -----

  clean_peptide <- function(x){
    str_remove(x, "\\..*") %>%
      str_remove_all("\\_") %>%
      str_remove_all("\\[.*?\\]") %>%
      str_trim()
  }

  # start execution ------

  # read spectronaut report

  suppressMessages(
    suppressWarnings(

      spctrnt_report <- read_tsv(
        path_to_file,
        na = c(
          "NA",
          "Filtered")
      )
    ))

  spct_df_min <- spctrnt_report %>%
    dplyr::select(
      PG.UniProtIds,
      PEP.IsProteotypic,
      EG.PrecursorId,
      EG.ModifiedSequence,
      contains(".d.EG.TotalQuantity")) %>%
    rename_with(~str_remove(.x, ".d.EG.TotalQuantity \\(Settings\\)")) %>%
    janitor::clean_names() %>%
    dplyr::rename(protein = pg_uni_prot_ids) %>%
    # define peptides that can appear more than once in a protein sequence
    mutate(
      multi_peptide = case_when(
        str_count(protein, "\\,") > 0 ~ TRUE,
        TRUE ~ FALSE
      ),
      # annotate N-termini (only acetyl supported for now)
      nterm_modif = case_when(
        str_detect(eg_modified_sequence, "Acetyl \\(N-term\\)") ~ "Acetyl",
        TRUE ~ "n"
        # no support for dimethyl yet, but can be implemented here needed
      )
    ) %>%
    # filter for proteotypic peptides
    filter(pep_is_proteotypic == proteotypic) %>%
    # exclude peptides with multiple appearances in a protein sequence
    filter(!multi_peptide) %>%
    # strip peptide sequence
    mutate(
      peptide = clean_peptide(eg_precursor_id)
    ) %>%
    mutate(
      nterm_modif_peptide = paste(
        nterm_modif,
        peptide,
        sep = "_"
      )) %>%
    dplyr::select(
      nterm_modif_peptide,
      nterm_modif,
      peptide,
      protein,
      where(is.double)
    )

  pept2protein <- spct_df_min %>%
    dplyr::select(
      nterm_modif_peptide,
      protein
    ) %>%
    distinct()

  # summ up the abundances for the same nterminally modified peptide
  spct_df_min_sum <- spct_df_min %>%
    group_by(nterm_modif_peptide, nterm_modif, peptide) %>%
    summarise(across(
      where(is.double),
      sum,
      na.rm = TRUE
    )) %>%
    ungroup() %>%
    mutate(
      across(
        where(is.double),
        ~log2(.x)
      )
    ) %>%
    left_join(
      pept2protein,
      by = "nterm_modif_peptide",
      multiple = "first"
    ) %>%
    dplyr::arrange(peptide) %>%
    dplyr::relocate(
      nterm_modif_peptide,
      nterm_modif,
      peptide,
      protein
    )

  return(spct_df_min_sum)
}
#' Process report file from DIANN search results for further analysis with TermineR
#'
#' @title diann_adapter
#'
#' @param path_to_file path to report file with DIANN search results in tsv format
#' @param proteotypic keep only proteotypic peptides. Default is TRUE.
#' @param summarization quantitative summarization approach for features identified by several precursors.
#' "SUM" will sum up the intensities of the precursors. "MAX" will keep the intensity of the precursor with higher intensity.
#' Default is "SUM".
#'
#' @format A data frame with at least 4 columns
#' \describe{
#'  \item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
#'  \item{nterm_modif}{N-terminal modification. Currently no N-terminal modification is supported for DIANN data}
#'  \item{peptide}{Peptide sequence}
#'  \item{protein}{Protein ID based on Uniprot Accession}
#' }
#'
#' @description
#' Process report file from DIANN search results for further analysis with TermineR
#'
#' @importFrom diann diann_load diann_matrix
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr select mutate distinct arrange group_by summarize ungroup
#' @importFrom tidyselect where
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
#' @importFrom data.table melt dcast
#'
#' @export
#' @author Miguel Cosenza-Contreras

diann_adapter <- function(
    path_to_file,
    proteotypic = TRUE,
    summarization = "SUM" # options: "SUM" or "MAX"
    ) {

  require(diann)
  require(tibble)
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(data.table)
  
  data("unimod_id_to_name_mapping", package = "TermineR")
  
  # the following function modifications are based on the diann R package code hosted on GitHub
  # https://github.com/vdemichev/diann-rpackage/blob/master/R/diann-R.R
  
  # pivot_aggregate function modified to sum up instead of get max value
  pivot_aggregate_2 <- function(
    df, 
    sample.header, 
    id.header, 
    quantity.header) {
    
    x <- melt.data.table(
      df, 
      id.vars = c(
        sample.header, 
        id.header), 
      measure.vars = c(
        quantity.header))
    
    x$value[which(x$value == 0)] <- NA
    
    piv <- as.data.frame(
      dcast.data.table(x, 
                       as.formula(paste0(id.header,'~',sample.header)), 
                       value.var = "value", 
                       fun.aggregate = function(x) sum(x, na.rm=TRUE))) 
    
    rownames(piv) <- piv[[1]]
    
    piv[[1]] <- NULL
    
    piv <- piv[order(rownames(piv)),]
    
    piv = as.matrix(piv)
    
    piv[is.infinite(piv)] <- NA
    
    # check if value is 0
    piv[piv == 0] <- NA
    
    piv
  }
  
  # modify diann_matrix function to sum up instead of calculating max
  diann_matrix_2 <- function(
    x, 
    id.header = "Precursor.Id", 
    quantity.header = "Precursor.Normalised", 
    proteotypic.only = F, 
    q = 0.01, 
    protein.q = 1.0, 
    pg.q = 1.0, 
    gg.q = 1.0) {
    
    df <- as.data.table(x)
    
    if (proteotypic.only) df <- df[which(df[['Proteotypic']] != 0),]
    
    df <- unique(df[which(df[[id.header]] != "" & df[[quantity.header]] > 0 & df[['Q.Value']] <= q & df[['Protein.Q.Value']] <= protein.q & df[['PG.Q.Value']] <= pg.q & df[['GG.Q.Value']] <= gg.q),c("File.Name", id.header, quantity.header),with=FALSE])
    
    is_duplicated = any(duplicated(paste0(df[["File.Name"]],":",df[[id.header]])))
    
    if (is_duplicated) {
      
      warning("Multiple quantities per id: the sum of these will be calculated")
      pivot_aggregate_2(df,"File.Name",id.header,quantity.header)
      
    } else {
      
      pivot(df,"File.Name",id.header,quantity.header)
      
    }
  }
  

  diann_df_min <- diann_load(path_to_file) %>%
    # extract modifications from modified_sequence check for everything inside '()'
    mutate(
      first_modif = str_extract(Modified.Sequence, "\\((.*?)\\)") %>% 
                                map_chr(~paste(.x, collapse = ";"))
    ) %>%
    # substitute the modification, including the brackets, with a Z
    mutate(
      min_first_mod_seq = str_replace(Modified.Sequence, "\\((.*?)\\)", "Z")
    ) %>%
    # get the location of those modifications
    mutate(
      # get the start of the first modification only using str_locate
      first_modif_locat = str_locate(min_first_mod_seq, "Z")[, "start"],
      id_nr = parse_number(first_modif)) %>%
    left_join(.,
              unimod_id_to_name_mapping) %>%
    mutate(
      nterm_modif = case_when(
        is.na(name) ~ "n",
        first_modif_locat == 1 & !is.na(name) ~ name,
        first_modif_locat != 1 & !is.na(name) ~ "n"
      )
    ) %>%
    mutate(
      nterm_modif_peptide = paste(
        nterm_modif,
        Stripped.Sequence,
        sep = "_"),
    ) %>%
    relocate(
      nterm_modif_peptide,
      nterm_modif,
      .before = Modified.Sequence
    )
  
  peptide2protein <- diann_df_min %>%
    dplyr::select(
      protein = Protein.Ids,
      peptide = Stripped.Sequence,
      nterm_modif_peptide,
      nterm_modif,
    ) %>%
    mutate(
      protein = str_remove(
        protein,
        " .*")
    ) %>%
    dplyr::select(
      protein,
      peptide,
      nterm_modif_peptide
    ) %>%
    distinct()
  
  # evaluate if summarization == "SUM" or "MAX" and choose the right diann_amtrix function
  
  if(summarization == "SUM"){
  
    diann_df_mat_min <- diann_matrix_2(
      diann_df_min,
      id.header = "nterm_modif_peptide",
      quantity.header = "Precursor.Normalised",
      proteotypic.only = TRUE,
      q = 0.01)
  
  } else if(summarization == "MAX"){
  
    diann_df_mat_min <- diann_matrix(
      diann_df_min,
      id.header = "nterm_modif_peptide",
      quantity.header = "Precursor.Normalised",
      proteotypic.only = TRUE,
      q = 0.01)
  
  }

  colnames(diann_df_mat_min) <- str_remove(
    colnames(diann_df_mat_min),
    ".*\\\\") %>%
    str_remove(
      "\\..*")

  diann_df_mat_min <- diann_df_mat_min %>%
    as.data.frame() %>%
    rownames_to_column("nterm_modif_peptide") %>%
    arrange(nterm_modif_peptide) %>%
    # apply log2 transformation
    mutate(
      across(
        where(is.double),
        ~log2(.x)
      )
    )
      
  # apply MAD scaling
  scaled_diann_df_mat_min <- diann_df_mat_min %>%
    # transform to long format to facilitate calculations at different levels
    pivot_longer(cols = where(is.double),
                 names_to = "sample",
                 values_to = "log2_ratio_2_ref") %>%
    
    # calculate median ratios per desired level, per sample
    group_by(nterm_modif_peptide, sample) %>%
    summarize(log2_rat2ref_group = median(log2_ratio_2_ref,
                                          na.rm = TRUE),
              .groups = "drop") %>%
    # calculate median ratios per sample
    group_by(sample) %>%
    mutate(Mi = median(log2_rat2ref_group,
                       na.rm = TRUE)) %>%
    ungroup() %>%
    
    # global median across all samples
    mutate(M0 = median(Mi,
                       na.rm = TRUE)) %>%
    
    # median centered ratios based on median ratios per sample
    mutate(RCij = log2_rat2ref_group - Mi) %>%
    
    # calculate median absolute deviation per sample
    group_by(sample) %>%
    mutate(MADi = median(abs(RCij),
                         na.rm = TRUE)) %>%
    ungroup() %>%
    
    # calculate global median absolute deviation
    mutate(MAD0 = median(MADi,
                         na.rm = TRUE)) %>%
    ungroup() %>%
    
    # calculate scaled ratios
    mutate(RNij = (RCij / MADi) * MAD0 + M0)
  
  # turn back to wide format using RNij as abundance
  scaled_diann_df_mat_min <- scaled_diann_df_mat_min %>%
    dplyr::select(
      nterm_modif_peptide,
      sample,
      RNij
    ) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "RNij"
    ) %>%
    separate(
      col = nterm_modif_peptide,
      into = c("nterm_modif",
               "peptide"),
      sep = "_",
      remove = FALSE
    ) %>%
    left_join(
      peptide2protein,
      multiple = "first"
    ) %>%
    relocate(
      nterm_modif_peptide,
      nterm_modif,
      peptide,
      protein
    ) %>%
    arrange(peptide)

  return(scaled_diann_df_mat_min)
}
#' Process psm.tsv files from FragPipe label-free search results for further analysis with TermineR
#'
#' @title fragpipe_lf_adapter
#' @param parent_dir directory of FragPipe search results
#' @param annotation_file_path path to the annotation file with sample names. It should contain at least two colums: `run`, matching the name of the LC-MS/MS run, and `sample_name`, with a readable sample identifier.
#' @param grouping_var define the grouping variable for MAD scaling. Default is "nterm_modif_peptide" for normalization at the feature level of "Peptide + N-terminal modification".
#'
#' @format A data frame with at least 4 columns
#' \describe{
#'  \item{nterm_modif_peptide}{Peptide identification merging N-terminal modification + peptide sequence}
#'  \item{nterm_modif}{N-terminal modification}
#'  \item{peptide}{Peptide sequence}
#'  \item{protein}{Protein ID based on Uniprot Accession}
#' }
#'
#'
#' @description
#' Process psm.tsv files from FragPipe search results for further analysis with TermineR
#'
#' @importFrom dplyr mutate filter select bind_rows left_join group_by summarize ungroup arrange distinct slice_max
#' @importFrom tidyr pivot_wider separate
#' @importFrom stringr str_detect str_remove
#' @importFrom readr read_tsv
#' @importFrom here here
#' @importFrom magrittr %>%
#'
#' @export
#' @author Miguel Cosenza-Contreras
fragpipe_lf_adapter <- function(
  parent_dir,
  annotation_file_path,
  grouping_var = "nterm_modif_peptide") {

  interest_cols <- c("Spectrum File",
                     "Peptide",
                     "Modified Peptide",
                     "PeptideProphet Probability",
                     "Intensity",
                     "Assigned Modifications",
                     "Is Unique",
                     "Protein ID",
                     "Gene")

# Construct the path to the psm.tsv files

if(!any(str_detect(list.files(parent_dir), "psm.tsv"))){

  # get the folder names of main project diretory
  folders_dir <- list.dirs(parent_dir,
                           full.names = TRUE,
                           recursive = FALSE)
  intern_foldrs <- list.dirs(parent_dir,
                             full.names = FALSE,
                             recursive = FALSE)
  psm_file_path <- paste0(folders_dir, "/psm.tsv")

} else {

  intern_foldrs <- "exp_1"

  folders_dir <- parent_dir

  psm_file_path <- paste0(parent_dir, "/psm.tsv")

}

# get a list of psm.tsv files
message("Loading psm.tsv files...")
list_psms <- purrr::map(.x = psm_file_path,
                        .f = read_tsv,
                        .progress = TRUE) %>%
                        suppressMessages()

# define names of list of psm files
names(list_psms) <- intern_foldrs

# combine all psm files into one
psm_tsv <- bind_rows(list_psms, .id = "experiment")

annotation_txt <- read_tsv(annotation_file_path) %>%
  na.omit()

psm_tsv_sel <- psm_tsv %>%
    dplyr::select(all_of(interest_cols)) %>%
    mutate(run = basename(`Spectrum File`)) %>%
    mutate(run = str_remove(run, ".pep.xml")) %>%
    mutate(run = str_remove(run, "interact-")) %>%
    dplyr::select(-`Spectrum File`) %>%
    left_join(annotation_txt) %>%
    mutate(
      nterm_modif = case_when(
        str_detect(`Assigned Modifications`, "N-term\\(304.207[0-9]\\)") ~ "TMT",
        str_detect(`Assigned Modifications`, "N-term\\(229.162[0-9]\\)") ~ "TMT",
        str_detect(`Assigned Modifications`, "N-term\\(42.010[0-9]\\)") ~ "Acetyl",
        str_detect(`Assigned Modifications`, "N-term\\(28.031[0-9]\\)") ~ "Dimethyl",
        str_detect(`Assigned Modifications`, "N-term\\(36.075[0-9]\\)") ~ "Dimethyl",
        str_detect(`Assigned Modifications`, "N-term\\(34.063[0-9]\\)") ~ "Dimethyl",
        TRUE ~ "n"
      )
    ) %>%
    mutate(
      nterm_modif_peptide = paste(
        nterm_modif,
        Peptide,
        sep = "_"
      )
    ) %>%
    dplyr::rename(
      peptide = Peptide,
      protein = `Protein ID`,
      intensity = Intensity
    ) %>%
    dplyr::select(
      nterm_modif_peptide,
      nterm_modif,
      intensity,
      peptide,
      protein,
      sample)

  scaled_ratios <- psm_tsv_sel %>%
    mutate(log2_intensity = log2(intensity)) %>%
    group_by(.data[[grouping_var]], sample) %>%
    summarize(log2_rat2ref_group = median(log2_intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(Mi = median(log2_rat2ref_group, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(M0 = median(Mi, na.rm = TRUE)) %>%
    mutate(RCij = log2_rat2ref_group - Mi) %>%
    group_by(sample) %>%
    mutate(MADi = median(abs(RCij), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(MAD0 = median(MADi, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(RNij = (RCij / MADi) * MAD0 + M0) %>%
    mutate(RNij = case_when(
      RNij == -Inf ~ NA,
      RNij == Inf ~ NA,
      TRUE ~ RNij
    ))

  pept2prot2modif <- psm_tsv_sel %>%
    dplyr::select(
      nterm_modif_peptide,
      nterm_modif,
      peptide,
      protein) %>%
    distinct()

  df_final_abundance <- scaled_ratios %>%
    dplyr::select(all_of(c(grouping_var, "sample", "RNij"))) %>%
    mutate(
      sample = str_remove(sample, "ratio_"),
    ) %>%
    pivot_wider(
      names_from = "sample",
      values_from = "RNij"
    ) %>%
    left_join(pept2prot2modif, by = grouping_var) %>%
    dplyr::relocate(
      nterm_modif_peptide,
      peptide,
      nterm_modif,
      protein
    )

  return(df_final_abundance)
}
