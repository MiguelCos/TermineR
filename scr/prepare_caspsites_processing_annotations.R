# LOAD REQUIRED PACKAGES -----

library(here)
library(tidyverse)

# PREPARE CASPSITES PROCESSING ANNOTATIONS -----

caspsites_search_url <- "https://www.caspsites.org/search/output"
caspsites_download_url <- "https://www.caspsites.org/search/download-sites"

caspsites_datasets <- c(
  "Caspase-1",
  "Caspase-2",
  "Caspase-3",
  "Caspase-4",
  "Caspase-5",
  "Caspase-6",
  "Caspase-7",
  "Caspase-8",
  "Caspase-9",
  "DegraBase",
  "1541",
  "Fas",
  "JG98",
  "Untreated Lysate"
)

caspsites_raw_dir <- here("data-raw/caspsites")
dir.create(caspsites_raw_dir, recursive = TRUE, showWarnings = FALSE)

caspsites_raw_csv <- here(caspsites_raw_dir, "caspsites_all_sites.csv")
caspsites_processing_tsv <- here(caspsites_raw_dir, "human_caspsites_processing.tsv")

encode_form <- function(form_names, form_values){

  paste(
    utils::URLencode(form_names, reserved = TRUE),
    utils::URLencode(form_values, reserved = TRUE),
    sep = "=",
    collapse = "&"
  )

}

download_caspsites_csv <- function(output_file){

  form_names <- c(
    rep("protease-list", length(caspsites_datasets)),
    "search-protein",
    "cleavageSite",
    "p4-residue",
    "p3-residue",
    "p2-residue",
    "p1-residue",
    "p1'-residue",
    "p2'-residue",
    "p3'-residue",
    "p4'-residue"
  )

  form_values <- c(
    caspsites_datasets,
    "",
    "custom-residues",
    rep("X", 8)
  )

  post_body <- encode_form(form_names, form_values)

  search_html <- tempfile("caspsites_search", fileext = ".html")
  handle <- httr::handle("https://www.caspsites.org")

  search_response <- try(
    httr::POST(
      url = caspsites_search_url,
      body = post_body,
      encode = "raw",
      httr::content_type("application/x-www-form-urlencoded"),
      httr::write_disk(search_html, overwrite = TRUE),
      httr::timeout(300),
      httr::config(connecttimeout = 60),
      handle = handle
    ),
    silent = TRUE
  )

  if(inherits(search_response, "try-error")){

    warning(
      "CaspSites search response ended before R finished reading the HTML page; ",
      "continuing if the CSV download endpoint is still populated."
    )

  } else if(httr::http_error(search_response)){

    stop(
      "CaspSites search request failed with HTTP status ",
      httr::status_code(search_response),
      "."
    )

  }

  if(!file.exists(search_html) || file.info(search_html)$size == 0){

    stop(
      "CaspSites search request did not return a usable response."
    )

  }

  download_response <- httr::GET(
    url = caspsites_download_url,
    httr::write_disk(output_file, overwrite = TRUE),
    httr::timeout(300),
    httr::config(connecttimeout = 60),
    handle = handle
  )

  if(httr::http_error(download_response)){

    stop(
      "CaspSites CSV download failed with HTTP status ",
      httr::status_code(download_response),
      "."
    )

  }
  invisible(output_file)

}

collapse_values <- function(values){

  values <- values[!is.na(values)]
  values <- stringr::str_squish(as.character(values))
  values <- values[values != "" & values != "None"]

  if(length(values) == 0){

    return("")

  }

  paste(sort(unique(values)), collapse = "|")

}

download_caspsites_csv(caspsites_raw_csv)

caspsites_raw <- readr::read_csv(
  caspsites_raw_csv,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = FALSE
)

human_caspsites_processing <- caspsites_raw %>%
  transmute(
    dataset = .data[["Dataset"]],
    protein = .data[["Accession"]],
    caspsites_uniprot_entry_name = .data[["Uniprot ID"]],
    caspsites_p1_prime_position = suppressWarnings(as.integer(.data[["Cleavage Site"]])),
    caspsites_p4 = .data[["P4"]],
    caspsites_p3 = .data[["P3"]],
    caspsites_p2 = .data[["P2"]],
    caspsites_p1 = .data[["P1"]],
    caspsites_p1_prime = .data[["P1'"]],
    caspsites_p2_prime = .data[["P2'"]],
    caspsites_p3_prime = .data[["P3'"]],
    caspsites_p4_prime = .data[["P4'"]],
    caspsites_peptide = .data[["Peptide"]],
    caspsites_cell_line = .data[["Cell Line"]],
    caspsites_perturbagen = .data[["Perturbagen"]],
    caspsites_label = .data[["Label"]],
    caspsites_pubmed = .data[["Pubmed"]],
    caspsites_doi = .data[["DOI"]],
    caspsites_author = .data[["Author"]]
  ) %>%
  filter(
    !is.na(protein),
    !is.na(caspsites_p1_prime_position)
  ) %>%
  mutate(
    caspsites_p1_position = caspsites_p1_prime_position - 1L,
    caspsites_cleavage_sequence = paste0(
      caspsites_p4,
      caspsites_p3,
      caspsites_p2,
      caspsites_p1,
      caspsites_p1_prime,
      caspsites_p2_prime,
      caspsites_p3_prime,
      caspsites_p4_prime
    ),
    caspsites_cleavage_site = paste0(
      caspsites_p4,
      caspsites_p3,
      caspsites_p2,
      caspsites_p1,
      "|",
      caspsites_p1_prime,
      caspsites_p2_prime,
      caspsites_p3_prime,
      caspsites_p4_prime
    )
  ) %>%
  group_by(
    protein,
    caspsites_p1_prime_position
  ) %>%
  summarise(
    caspsites_p1_position = dplyr::first(caspsites_p1_position),
    caspsites_cleavage_sequence = collapse_values(caspsites_cleavage_sequence),
    caspsites_cleavage_site = collapse_values(caspsites_cleavage_site),
    caspsites_datasets = collapse_values(dataset),
    caspsites_uniprot_entry_names = collapse_values(caspsites_uniprot_entry_name),
    caspsites_peptides = collapse_values(caspsites_peptide),
    caspsites_cell_lines = collapse_values(caspsites_cell_line),
    caspsites_perturbagens = collapse_values(caspsites_perturbagen),
    caspsites_labels = collapse_values(caspsites_label),
    caspsites_pubmed = collapse_values(caspsites_pubmed),
    caspsites_doi = collapse_values(caspsites_doi),
    caspsites_authors = collapse_values(caspsites_author),
    caspsites_evidence_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(
    protein,
    caspsites_p1_prime_position
  )

readr::write_tsv(
  human_caspsites_processing,
  caspsites_processing_tsv
)

save(
  human_caspsites_processing,
  file = here("data/human_caspsites_processing.rda"),
  compress = "xz",
  version = 2
)
