annotate_canonical_start_fixture <- function(
  protein_sequences,
  peptides_df,
  sense,
  specificity
) {
  fasta_path <- tempfile(
    fileext = ".fasta"
  )

  on.exit(
    unlink(fasta_path)
  )

  fasta_headers <- paste0(
    ">sp|",
    names(protein_sequences),
    "|CANONICAL_START_FIXTURE"
  )

  writeLines(
    as.vector(rbind(fasta_headers, protein_sequences)),
    con = fasta_path
  )

  peptides_df$nterm_modif_peptide <- paste0(
    peptides_df$nterm_modif,
    "_",
    peptides_df$peptide
  )

  peptides_df$sample1 <- 1

  annotated <- annotate_neo_termini(
    peptides_df = peptides_df,
    fasta_location = fasta_path,
    sense = sense,
    specificity = specificity,
    organism = "human",
    pssm_prediction = FALSE
  )

  expected_feature_keys <- paste(
    peptides_df$protein,
    peptides_df$nterm_modif_peptide
  )

  annotated_feature_keys <- paste(
    annotated$protein,
    annotated$nterm_modif_peptide
  )

  expect_equal(
    nrow(annotated),
    nrow(peptides_df)
  )

  expect_setequal(
    annotated_feature_keys,
    expected_feature_keys
  )

  annotated[
    match(expected_feature_keys, annotated_feature_keys),
    ,
    drop = FALSE
  ]
}

test_that(
  "mapped acetylated and free starts take precedence over processing annotations",
  {
    # Packaged annotations place UniProt INIT_MET on A1L3X0 at residue 1
    # and a TargetP signal cleavage on A0A140TA83 at P1 position 3.
    protein_sequences <- c(
      A1L3X0 = "MACDEKFGHILMNPQ",
      A0A140TA83 = "DQFGHATCIRKMN",
      TESTFREE = "PELTQKAYHVGCR",
      TESTINTERNAL = "VVDWGVIKGLRPLTQ",
      TESTUNMAPPED = "WQLAPDKSFVREHT"
    )

    peptides_df <- data.frame(
      protein = c(
        "A1L3X0",
        "A1L3X0",
        "A0A140TA83",
        "TESTFREE",
        "A1L3X0",
        "A1L3X0",
        "TESTINTERNAL",
        "TESTUNMAPPED"
      ),
      peptide = c(
        "MACDEK",
        "MACDEK",
        "DQFGHA",
        "PELTQK",
        "ACDEK",
        "MACDE",
        "DWGVIK",
        "MMMMK"
      ),
      nterm_modif = c("Acetyl", "n", "n", "n", "n", "TMT", "Acetyl", "n")
    )

    annotated <- annotate_canonical_start_fixture(
      protein_sequences = protein_sequences,
      peptides_df = peptides_df,
      sense = "C",
      specificity = "K|R"
    )

    expect_equal(
      annotated$peptide_start,
      c(1, 1, 1, 1, 2, 1, 3, NA)
    )

    expect_identical(
      annotated$nterm_modif,
      peptides_df$nterm_modif
    )

    expect_identical(
      annotated$specificity[1:4],
      c("specific", "specific", "semi_Nterm", "specific")
    )

    expect_identical(
      annotated$uniprot_processing_type[c(1, 2, 5, 6)],
      rep("INIT_MET", 4)
    )

    expect_identical(
      annotated$processing_type_targetp[3],
      "SIGNAL"
    )

    expect_identical(
      annotated$uniprot_processing_type[4],
      "not_canonical_no_procc_annot"
    )

    expect_identical(
      annotated$processing_type_terminer,
      c(
        "Intact_ORF_Acetylated",
        "Intact_ORF_free",
        "Intact_ORF_free",
        "Intact_ORF_free",
        "INIT_MET",
        "INIT_MET",
        "not_canonical",
        "not_canonical_no_procc_annot"
      )
    )

    expect_identical(
      as.character(annotated$processing_type),
      c("P1", "P1", "P1", "P1", "P2", "P2", "not_canonical", "not_canonical")
    )

    expect_identical(
      levels(annotated$processing_type),
      c("P1", "P2", "SP", "MTS", "cTP", "luTP", "TRANSIT", "PROPEP", "not_canonical")
    )

    expect_identical(
      annotated$source_processing_annotation,
      c("Other", "Other", "Other", "Other", "Uniprot", "Uniprot", "Other", "Other")
    )
  }
)

test_that(
  "an exact peptide start defines P1 even when specificity selects the C-terminal bond",
  {
    peptides_df <- data.frame(
      protein = "TESTSEMIC",
      peptide = "DQLMAG",
      nterm_modif = "Acetyl"
    )

    annotated <- annotate_canonical_start_fixture(
      protein_sequences = c(TESTSEMIC = "DQLMAGVLKPRFTM"),
      peptides_df = peptides_df,
      sense = "N",
      specificity = "D|E"
    )

    expect_equal(
      unname(annotated$peptide_start),
      1L
    )

    expect_equal(
      unname(annotated$p1_prime_position),
      7L
    )

    expect_identical(
      annotated$specificity,
      "semi_Cterm"
    )

    expect_identical(
      annotated$processing_type_terminer,
      "Intact_ORF_Acetylated"
    )

    expect_identical(
      as.character(annotated$processing_type),
      "P1"
    )
  }
)
