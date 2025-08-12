# LOAD REQUIRED PACKAGES -----

library(here)
library(tidyverse)
library(drawProteins)
library(seqinr)

# Uniprot processing features ----

mol_processing_feat <- c("INIT_MET",
                         "PROPEP",
                         "SIGNAL",
                         "TRANSIT",
                         "CHAIN")

# PREPARE HUMAN ANNOTATION ---------------------------------------------------

## fasta location

fasta_location <- here("data-raw/uniprotkb_Human_AND_model_organism_9606_2024_05_07.fasta")


## read and load fasta

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

## get features
if(!file.exists(here("data-raw/human_uniprot_features_pkd.rds"))){

  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)

  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
   
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/human_uniprot_features_pkd.rds"))
    
  
} else {
  
  uniprot_features <- read_rds(here("data-raw/human_uniprot_features_pkd.rds"))

}

## convert to dataframe

df_uniprot_human_features <- purrr::map(successful_results,
                                  drawProteins::feature_to_dataframe)

df_human_features <- bind_rows(df_uniprot_human_features) #%>%

write_tsv(df_human_features, here("data-raw/human_uniprot_processing_features.tsv"))

# PREPARE HUMAN + ISOFORMS ANNOTATION ---------------------------------------------------

## fasta location

fasta_location <- here("data-raw/prep_annotations/human_uniprot_swissprot_canon_n_isoforms_2024_08.fasta")


## read and load fasta

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

## get features
if(!file.exists(here("data-raw/prep_annotations/human_isoform_uniprot_features_pkd.rds"))){
  
  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)
  
  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
  
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/prep_annotations/human_isoform_uniprot_features_pkd.rds"))
  
  
} else {
  
  successful_results <- read_rds(here("data-raw/prep_annotations/human_isoform_uniprot_features_pkd.rds"))
  
}

## convert to dataframe

safe_feature_to_dataframe <- purrr::safely(drawProteins::feature_to_dataframe)

df_uniprot_human_iso_features <- purrr::map(successful_results,
                                            safe_feature_to_dataframe)

# check each df_uniprot_human_iso_features[[i]]$result for null values

df_uniprot_human_iso_features <- df_uniprot_human_iso_features %>%
  purrr::keep(~!is.null(.x$result))

# eliminate the 'error' sub-element from each element of the list

df_uniprot_human_iso_features <- purrr::map(df_uniprot_human_iso_features, ~.x$result)

df_human_iso_features <- bind_rows(df_uniprot_human_iso_features) #%>%
  # change colnames; eliminate result$ from the start of each colname
  rename_all(~str_remove(.x, "result\\$"))

write_tsv(
  df_human_iso_features, 
  here("data-raw/prep_annotations/human_iso_uniprot_processing_features.tsv"))

# PREPARE PIG ANNOTATION ---------------------------------------------------

## fasta location

fasta_location <- here("data-raw/prep_annotations/sus_scrofa_pig_UP000008227_9823.fasta")

## read and load fasta

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

## get features

if(!file.exists(here("data-raw/pig_uniprot_features_pkd.rds"))){

  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)

  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
   
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/prep_annotations/pig_uniprot_features_pkd.rds"))
    
  
} else {
  
  uniprot_features <- read_rds(here("data-raw/prep_annotations/pig_uniprot_features_pkd.rds"))

}

## convert to dataframe

df_uniprot_pig_features <- purrr::map(successful_results,
                                  drawProteins::feature_to_dataframe)

df_pig_features <- bind_rows(df_uniprot_pig_features) #%>%

write_tsv(df_pig_features, here("data-raw/prep_annotations/pig_uniprot_processing_features.tsv"))

# PREPARE ARABIDOPSIS THALIANA ANNOTATION -----

## read and load fasta ----

## fasta location

fasta_location <- here("data-raw/prep_annotations/uniprotkb_Arabidopsis_thaliana_AND_mode_2024_09_12_fas.fasta")

## read and load fasta

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

## get features ----

if(!file.exists(here("data-raw/prep_annotations/arabidopsis_uniprot_features_pkd.rds"))){

  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)

  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
   
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/prep_annotations/arabidopsis_uniprot_features_pkd.rds"))
    
  
} else {
  
  successful_results <- read_rds(here("data-raw/prep_annotations/arabidopsis_uniprot_features_pkd.rds"))

}

## convert to dataframe ---- 
safe_feature_to_dataframe <- purrr::safely(drawProteins::feature_to_dataframe)

df_uniprot_arabidopsis_features <- purrr::map(successful_results,
                                              safe_feature_to_dataframe)

# extract the $result subelement of each element of the list

df_uniprot_arabidopsis_features <- purrr::map(df_uniprot_arabidopsis_features, ~.x$result)

# exclude any NULL element of the list

df_uniprot_arabidopsis_features <- df_uniprot_arabidopsis_features %>%
  purrr::keep(~!is.null(.x))

df_arabidopsis_features <- bind_rows(df_uniprot_arabidopsis_features) %>%
  filter(type %in% mol_processing_feat)
  
write_tsv(df_arabidopsis_features, 
          here("data-raw/prep_annotations/arabidopsis_uniprot_processing_features.tsv"))

# PREPARE E. COLI UNIPROT ANNOTATION ----- 

## read and load fasta ----

## fasta location

fasta_location <- here("data-raw/prep_annotations/uniprotkb_e_coli_AND_model_organism_833_2024_09_13.fasta")

## read and load fasta

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

## get features ----

if(!file.exists(here("data-raw/prep_annotations/e_coli_uniprot_features_pkd.rds"))){

  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)

  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
   
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/prep_annotations/e_coli_uniprot_features_pkd.rds"))
    
  
} else {
  
  successful_results <- read_rds(here("data-raw/prep_annotations/e_coli_uniprot_features_pkd.rds"))

}

## convert to dataframe ----

safe_feature_to_dataframe <- purrr::safely(drawProteins::feature_to_dataframe)

df_uniprot_e_coli_features <- purrr::map(successful_results,
                                              safe_feature_to_dataframe)

# extract the $result subelement of each element of the list

df_uniprot_e_coli_features <- purrr::map(df_uniprot_e_coli_features, ~.x$result)

# exclude any NULL element of the list

df_uniprot_e_coli_features <- df_uniprot_e_coli_features %>%
  purrr::keep(~!is.null(.x))

df_e_coli_features <- bind_rows(df_uniprot_e_coli_features) %>%
  filter(type %in% mol_processing_feat)

write_tsv(df_e_coli_features, 
          here("data-raw/prep_annotations/e_coli_uniprot_processing_features.tsv"))


# PREPARE C. ELEGANS UNIPROT ANNOTATION ----- 

## read and load fasta ----

## fasta location

fasta_location <- here("data-raw/prep_annotations/uniprotkb_C_elegans_AND_model_organism_2025_08_12.fasta")

## read and load fasta

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

## get features ----

if(!file.exists(here("data-raw/prep_annotations/c_elegans_uniprot_features_pkd.rds"))){

  time1 <- Sys.time()
  safe_get_features <- purrr::safely(drawProteins::get_features)

  uniprot_features <- purrr::map(.x = fasta_df$protein,
                                 .f = safe_get_features,
                                 .progress = TRUE)
   
  # To get only the successful results and ignore the errors
  successful_results <- purrr::map(uniprot_features, "result")
  time2 <- Sys.time()
  
  write_rds(successful_results, 
            file = here("data-raw/prep_annotations/c_elegans_uniprot_features_pkd.rds"))
    
  
} else {
  
  successful_results <- read_rds(here("data-raw/prep_annotations/c_elegans_uniprot_features_pkd.rds"))

}

## convert to dataframe ----

safe_feature_to_dataframe <- purrr::safely(drawProteins::feature_to_dataframe)

df_uniprot_e_coli_features <- purrr::map(successful_results,
                                              safe_feature_to_dataframe)

# extract the $result subelement of each element of the list

df_uniprot_e_coli_features <- purrr::map(df_uniprot_e_coli_features, ~.x$result)

# exclude any NULL element of the list

df_uniprot_e_coli_features <- df_uniprot_e_coli_features %>%
  purrr::keep(~!is.null(.x))

df_e_coli_features <- bind_rows(df_uniprot_e_coli_features) %>%
  filter(type %in% mol_processing_feat)

write_tsv(df_e_coli_features, 
          here("data-raw/prep_annotations/c_elegans_uniprot_processing_features.tsv"))

