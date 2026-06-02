# LOAD REQUIRED PACKAGES -----

library(here)
library(tidyverse)

# PREPARE UNIMOD ID TO NAME MAPPING -----

unimod_obo_url <- "https://www.unimod.org/obo/unimod.obo"
unimod_obo_location <- here("data-raw/prep_annotations/unimod.obo")

dir.create(
  dirname(unimod_obo_location),
  recursive = TRUE,
  showWarnings = FALSE
)

download.file(
  url = unimod_obo_url,
  destfile = unimod_obo_location,
  mode = "wb"
)

parse_unimod_xref <- function(term_lines, xref_name){

  xref_pattern <- paste0("^xref: ", xref_name, " \"([^\"]*)\"")
  xref_line <- str_subset(term_lines, xref_pattern)

  if(length(xref_line) == 0){

    return(NA_character_)

  }

  str_match(xref_line[[1]], xref_pattern)[[2]]

}

parse_unimod_term <- function(term_lines){

  id <- str_subset(term_lines, "^id: UNIMOD:[0-9]+$") %>%
    str_remove("^id: ")

  name <- str_subset(term_lines, "^name: ") %>%
    str_remove("^name: ")

  if(length(id) == 0 || length(name) == 0){

    return(NULL)

  }

  tibble(
    id = id[[1]],
    id_nr = as.numeric(str_remove(id[[1]], "UNIMOD:")),
    name = name[[1]],
    monoisotopic_mass = as.numeric(parse_unimod_xref(term_lines, "delta_mono_mass")),
    average_mass = as.numeric(parse_unimod_xref(term_lines, "delta_avge_mass")),
    composition = parse_unimod_xref(term_lines, "delta_composition")
  )

}

unimod_obo_lines <- readLines(unimod_obo_location, warn = FALSE)

term_starts <- which(unimod_obo_lines == "[Term]")
term_ends <- c(term_starts[-1] - 1, length(unimod_obo_lines))

unimod_id_to_name_mapping <- map2(
  term_starts,
  term_ends,
  ~ parse_unimod_term(unimod_obo_lines[(.x + 1):.y])
) %>%
  compact() %>%
  bind_rows() %>%
  filter(
    id_nr > 0
  ) %>%
  arrange(
    id_nr
  )

write_tsv(
  unimod_id_to_name_mapping,
  here("data-raw/prep_annotations/unimod_id_to_name_mapping.tsv")
)

save(
  unimod_id_to_name_mapping,
  file = here("data/unimod_id_to_name_mapping.rda"),
  compress = "xz",
  version = 2
)
