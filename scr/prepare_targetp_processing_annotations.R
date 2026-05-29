# LOAD REQUIRED PACKAGES -----

library(here)
library(tidyverse)
library(janitor)

# PREPARE TARGETP PROCESSING ANNOTATIONS -----

process_targetp_summary <- function(targetp_location){

  read_tsv(
    targetp_location,
    skip = 1,
    show_col_types = FALSE
  ) %>%
    clean_names() %>%
    separate(
      number_id,
      into = c("database", "protein", "uniprot_entry_name"),
      sep = "\\|",
      remove = FALSE,
      fill = "right",
      extra = "merge"
    ) %>%
    dplyr::rename(
      protein_header = number_id
    ) %>%
    mutate(
      targetp_p1_position = str_remove(
        cs_position,
        "CS pos:\\s*"
      )
    ) %>%
    mutate(
      targetp_p1_position = str_remove(
        targetp_p1_position,
        "-.*"
      )
    ) %>%
    mutate(
      targetp_p1_position = as.numeric(targetp_p1_position)
    ) %>%
    filter(
      !is.na(targetp_p1_position)
    ) %>%
    dplyr::select(
      protein,
      targetp_category = prediction,
      targetp_p1_position
    ) %>%
    distinct(
      protein,
      .keep_all = TRUE
    )

}

targetp_files <- tribble(
  ~organism, ~targetp_location,
  "human", here("data-raw/targetp/human_uniprot_2024_09/h_sapiens_long_summary.targetp2"),
  "mouse", here("data-raw/targetp/mouse_uniprot_2025_01/mouse_uniprot_summary.targetp2"),
  "arabidopsis", here("data-raw/targetp/arabidopsis_thaliana/arabidopsis_thaliana_summary_2024_11.targetp2"),
  "rat", here("data-raw/targetp/rat_uniprot_2025_11/rat_reference_summary.targetp2"),
  "yeast", here("data-raw/targetp/yeast_uniprot_2026_05/yeast_reference_summary.targetp2"),
  "medicago_trucantula", here("data-raw/targetp/medicago_trucantula/medicago_trucantula_reference_summary.targetp2"),
  "rhizobium_melitoli", here("data-raw/targetp/rhizobium_melitoli/rhizobium_melitoli_reference_summary.targetp2"),
  "pig", here("data-raw/targetp/pig/pig_reference_summary.targetp2"),
  "human_iso", here("data-raw/targetp/human_iso/human_iso_reference_summary.targetp2"),
  "ecoli", here("data-raw/targetp/ecoli/ecoli_reference_summary.targetp2"),
  "c_elegans", here("data-raw/targetp/c_elegans/c_elegans_reference_summary.targetp2"),
  "synechocystis", here("data-raw/targetp/synechocystis/synechocystis_reference_summary.targetp2")
)

targetp_files <- targetp_files %>%
  filter(
    file.exists(targetp_location)
  )

if(nrow(targetp_files) == 0){

  stop("No TargetP summary files found in data-raw/targetp.")

}

for(i in seq_len(nrow(targetp_files))){

  organism <- targetp_files$organism[[i]]
  targetp_location <- targetp_files$targetp_location[[i]]
  targetp_processing <- process_targetp_summary(targetp_location)
  targetp_object_name <- paste0(organism, "_targetp_processing")

  write_tsv(
    targetp_processing,
    here("data-raw/targetp", paste0(targetp_object_name, ".tsv"))
  )

  assign(
    targetp_object_name,
    targetp_processing
  )

  save(
    list = targetp_object_name,
    file = here("data", paste0(targetp_object_name, ".rda")),
    compress = "xz",
    version = 2
  )

}
