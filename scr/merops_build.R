# data-raw/merops_build.R
# Build: merops_sites (known cleavages + window8), merops_pssm (per-protease PSSM)

library(tidyverse)
library(usethis)

# ---------- 0) Config ----------
ftp_dir <- "ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release"
fn_substrates <- file.path(ftp_dir, "substrates.txt")
fn_cleavage   <- file.path(ftp_dir, "cleavage.txt")           # name may be substrate_cut_sites in some releases

# Local temp
td <- tempdir()
f_substrates <- file.path(td, "substrates.txt")
f_cleavage   <- file.path(td, "cleavage.txt")

# ---------- 1) Download (FTP, anonymous) ----------
curl::curl_download(fn_substrates, f_substrates, mode = "wb")
curl::curl_download(fn_cleavage,   f_cleavage,   mode = "wb")

# ---------- 2) Read ----------
substrates <- read.delim(f_substrates, header = FALSE, quote = "\"",
                         na.strings = "\\N", stringsAsFactors = FALSE)
colnames(substrates) <- c("uniprot_acc", "sequence", "merops_substrate_id",
                          "description", "taxid", "comments")

cleavage <- read.delim(f_cleavage, header = FALSE, quote = "\"",
                       na.strings = "\\N", stringsAsFactors = FALSE)
# Column names are stable-ish but MEROPS exports vary slightly across releases.
# This covers the common layout you pasted:
colnames(cleavage) <- c("protease_id", "substrate_key", "position",
                        "evidence_type", "physio_context",
                        "substrate_name", "substrate_desc", "extra_annot",
                        "region_tested", "notes", "reference", "method", "comments")
cleavage$position <- suppressWarnings(as.integer(cleavage$position))

# ---------- 3) Join cleavage -> sequence ----------
# Try MEROPS internal substrate id first; if NA, try UniProt accession
cleav_seq <- cleavage %>%
  left_join(substrates %>% select(merops_substrate_id, sequence, uniprot_acc),
            by = c("substrate_key" = "merops_substrate_id")) %>%
  left_join(substrates %>% select(uniprot_acc, sequence) %>% rename(sequence2 = sequence),
            by = c("substrate_key" = "uniprot_acc")) %>%
  mutate(sequence = coalesce(sequence, sequence2),
         uniprot_acc = coalesce(uniprot_acc, substrate_key)) %>%
  select(-sequence2) %>%
  filter(!is.na(sequence), !is.na(position))

# ---------- 4) Build P4…P4′ windows ----------
AA20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
extract_window8 <- function(seq, pos, flank = 4) {
  # MEROPS convention: cleavage is between P1 (pos-1) and P1' (pos)
  # pos is the index of P1' in the protein sequence (1-based)
  start <- max(1, pos - flank)
  end   <- min(nchar(seq), pos + (flank - 1))
  w <- substr(seq, start, end)
  # Pad with X to length 8 (if near ends)
  if (nchar(w) < 8) w <- stringr::str_pad(w, 8, side = "right", pad = "X")
  if (nchar(w) < 8) w <- stringr::str_pad(w, 8, side = "left",  pad = "X")
  # Uppercase, non-standard letters to X
  w <- toupper(w)
  w <- chartr("BJOUZ*", "XXXXX", w)
  w
}

merops_sites <- cleav_seq %>%
  mutate(window8 = purrr::map2_chr(sequence, position, extract_window8, flank = 4)) %>%
  filter(!grepl("^X+$", window8)) %>%
  transmute(
    protease_id,
    uniprot_acc,
    position,                 # P1' index in UniProt sequence
    window8,
    evidence_type = coalesce(evidence_type, NA_character_),
    physio_context = coalesce(physio_context, NA_character_),
    method = coalesce(method, NA_character_),
    reference = coalesce(reference, NA_character_)
  ) %>%
  distinct()

# ---------- 5) Build PSSMs per protease ----------
# global background from all windows (AA frequency across all positions)
bg_from_windows <- function(windows) {
  vec <- strsplit(paste(windows, collapse = ""), "")[[1]]
  vec <- vec[vec %in% AA20]
  tab <- table(factor(vec, levels = AA20))
  as.numeric(tab) / sum(tab)
}

global_bg <- bg_from_windows(merops_sites$window8)

# Optional BLOSUM62-based prior (softened)
softmax <- function(x, t = 1) { ex <- exp(x / t); ex / sum(ex) }
make_blosum_prior <- function() {
  # Ship a tiny BLOSUM62 matrix in the package or depend on Biostrings if you prefer.
  # Below uses Biostrings; replace with an internal copy if you want to avoid Suggests.
  m <- as.matrix(Biostrings::blosum62[AA20, AA20])
  pri <- t(apply(m, 1, softmax, t = 3))  # temperature 3 = smoother prior
  dimnames(pri) <- list(AA20, AA20)
  pri
}
PRI <- tryCatch(make_blosum_prior(), error = function(e) NULL)

count_pos_matrix <- function(windows) {
  mats <- lapply(1:8, function(i) table(factor(substr(windows, i, i), levels = AA20)))
  do.call(cbind, mats) # 20 x 8
}

make_pssm <- function(windows, bg = global_bg, pc = 1.0, use_blosum_prior = TRUE) {
  if (!length(windows)) return(NULL)
  counts <- count_pos_matrix(windows) # 20x8
  if (use_blosum_prior && !is.null(PRI)) {
    # position-wise data-adaptive prior by mixing BLOSUM62 rows weighted by observed column frequencies
    counts_pc <- vapply(1:8, function(j) {
      col <- counts[, j, drop = FALSE]
      tot <- sum(col)
      w   <- if (tot > 0) as.numeric(col) / tot else rep(1/20, 20)
      # prior for the column = PRI %*% w
      prior_col <- PRI %*% w
      as.numeric(col + pc * prior_col)
    }, numeric(20))
  } else {
    counts_pc <- counts + pc
  }
  probs <- sweep(counts_pc, 2, colSums(counts_pc), "/")     # 20 x 8
  odds  <- sweep(probs, 1, bg, "/")                         # divide by background
  log2(odds)                                                # 20 x 8
}

protease_windows <- merops_sites %>%
  group_by(protease_id) %>%
  summarise(windows = list(window8), n_sites = dplyr::n(), .groups = "drop")

merops_pssm <- protease_windows %>%
  mutate(
    pssm = purrr::map(windows, make_pssm, bg = global_bg, pc = 1.0, use_blosum_prior = TRUE),
    bg   = list(setNames(global_bg, AA20))
  ) %>%
  select(-windows) %>%
  filter(!purrr::map_lgl(pssm, is.null))

# ---------- 6) Save in package ----------
usethis::use_data(merops_sites, merops_pssm, overwrite = TRUE)
