
library(here)
library(tidyverse)
library(drawProteins)
library(seqinr)

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


