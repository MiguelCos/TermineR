annotate_cleavage_residue_fixture <- function(
  protein_sequences,
  peptides_df,
  sense,
  specificity
) {
  fasta_file <- tempfile(
    fileext = ".fasta"
  )

  on.exit(
    unlink(fasta_file)
  )

  writeLines(
    as.vector(
      rbind(
        paste0(">sp|", peptides_df$protein, "|CLEAVAGE_FIXTURE"),
        protein_sequences
      )
    ),
    con = fasta_file
  )

  annotate_neo_termini(
    peptides_df = peptides_df,
    fasta_location = fasta_file,
    sense = sense,
    specificity = specificity,
    organism = "human",
    pssm_prediction = FALSE
  )
}

test_that(
  "cleavage residues match the annotated bond for both protease senses",
  {
    # Synthetic FASTA places the issue #42 sequence context at P1 position 195.
    issue_sequence <- paste0(
      strrep("A", 191),
      "TLLKDAAKVCREFT"
    )

    protein_sequences <- c(
      issue_sequence,
      issue_sequence,
      "DAAKVCREFT",
      "ATLLKDAAKVCR",
      "ATLLQAAKVCREFT",
      "ATLLKDAAKVCAQFT"
    )

    peptides_df <- data.frame(
      nterm_modif = c("TMT", "n", rep("TMT", 4)),
      peptide = c(rep("DAAKVCR", 4), "AAKVCR", "DAAKVCA"),
      protein = paste0("TEST", seq_along(protein_sequences)),
      sample1 = 1,
      stringsAsFactors = FALSE
    )

    peptides_df$nterm_modif_peptide <- paste0(
      peptides_df$nterm_modif,
      "_",
      peptides_df$peptide
    )

    expected <- data.frame(
      protein = peptides_df$protein,
      nterm_modif = peptides_df$nterm_modif,
      specificity = c(rep("specific", 4), "semi_Nterm", "semi_Cterm"),
      p1_position = c(195, 195, 0, 5, 5, 12),
      p1_prime_position = c(196, 196, 1, 6, 6, 13),
      p1_residue = c("K", "K", "N-term", "K", "Q", "A"),
      p1_prime_residue = c("D", "D", "D", "D", "A", "Q"),
      cleavage_site = c(
        "TLLK | DAAK",
        "TLLK | DAAK",
        "XXXX | DAAK",
        "TLLK | DAAK",
        "TLLQ | AAKV",
        "KVCA | QFTX"
      ),
      cleavage_sequence = c(
        "TLLKDAAK",
        "TLLKDAAK",
        "XXXXDAAK",
        "TLLKDAAK",
        "TLLQAAKV",
        "KVCAQFTX"
      )
    )

    for (sense in c("N", "C")) {
      specificity <- if (sense == "N") "D|E" else "K|R"

      # Batch the fixtures to avoid repeatedly loading the reference datasets.
      annotated <- annotate_cleavage_residue_fixture(
        protein_sequences = protein_sequences,
        peptides_df = peptides_df,
        sense = sense,
        specificity = specificity
      )

      expect_equal(
        as.data.frame(annotated[, names(expected)]),
        expected,
        info = paste("sense =", sense)
      )
    }
  }
)
