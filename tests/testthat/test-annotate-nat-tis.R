test_that(
  "NAT motifs and potential alternative starts use the mapped peptide N-terminus",
  {
    protein_sequences <- c(
      NAT_A_BOUNDARY = "AAAMVTQLKDTGFR",
      NAT_A_INTERNAL = "AAAAMVTQLKDTGFR",
      NAT_B = "AAAAAQMDAQKDEWPR",
      NAT_CEF_LEUCINE = "MLQGPKQTVFR",
      NAT_CEF_ALANINE = "MAFQKDTQPR",
      NAT_CEF_SERINE = "MSGTLKDTQPR",
      NAT_D = "MSGRGPAKDTLQFR",
      NAT_H_ASPARTATE = "AAAAAMDDDIAVKQTGFR",
      NAT_H_GLUTAMATE = "EEEIYQKPNLFR",
      ALT_START_NO_MOTIF = "AAAAAQMPPLVKDTPR",
      MET_AT_CTERM_BOND = "AAAAAADEEEIQKMTPR",
      FREE_PEPTIDE = "AAAAAMCTQPKDTVFR",
      UNKNOWN_MODIFICATION = "AAAAAMGTVLKDQPSR",
      UNMAPPED_PEPTIDE = "WQLAPDKSFVREHT"
    )

    peptides_df <- data.frame(
      protein = names(protein_sequences),
      peptide = c(
        "VTQLK",
        "VTQLK",
        "MDAQK",
        "MLQGPK",
        "MAFQK",
        "MSGTLK",
        "SGRGPAK",
        "DDDIAVK",
        "EEEIYQK",
        "MPPLVK",
        "DEEEIQK",
        "CTQPK",
        "GTVLK",
        "QFWLK"
      ),
      nterm_modif = c(rep("Acetyl", 11), "n", NA_character_, "Acetyl"),
      sample1 = 1
    )

    peptides_df$nterm_modif_peptide <- paste0(
      peptides_df$nterm_modif,
      "_",
      peptides_df$peptide
    )

    fasta_path <- tempfile(
      fileext = ".fasta"
    )

    on.exit(
      unlink(fasta_path)
    )

    fasta_headers <- paste0(
      ">sp|",
      names(protein_sequences),
      "|NAT_TIS_FIXTURE"
    )

    writeLines(
      as.vector(rbind(fasta_headers, protein_sequences)),
      con = fasta_path
    )

    annotated <- annotate_neo_termini(
      peptides_df = peptides_df,
      fasta_location = fasta_path,
      sense = "N",
      specificity = "D|E",
      organism = "human",
      pssm_prediction = FALSE
    )

    expect_equal(
      nrow(annotated),
      nrow(peptides_df)
    )

    expect_setequal(
      annotated$protein,
      peptides_df$protein
    )

    annotated <- annotated[
      match(peptides_df$protein, annotated$protein),
      ,
      drop = FALSE
    ]

    expect_identical(
      annotated$nterm_modif_peptide,
      peptides_df$nterm_modif_peptide
    )

    expected_annotations <- data.frame(
      peptide_start = c(5, 6, 7, 1, 1, 1, 2, 7, 1, 7, 7, 7, 7, NA),
      nat_match = c(
        "NatA",
        "NatA",
        "NatB",
        "NatC/E/F",
        "NatC/E/F",
        "NatC/E/F",
        "NatD",
        "NatH",
        "NatH",
        NA_character_,
        NA_character_,
        NA_character_,
        NA_character_,
        NA_character_
      ),
      potential_alt_tis = c(
        FALSE,
        TRUE,
        TRUE,
        FALSE,
        FALSE,
        FALSE,
        FALSE,
        TRUE,
        FALSE,
        TRUE,
        FALSE,
        FALSE,
        FALSE,
        FALSE
      ),
      row.names = peptides_df$protein
    )

    observed_annotations <- as.data.frame(
      annotated[, names(expected_annotations), drop = FALSE]
    )

    rownames(observed_annotations) <- annotated$protein

    expect_equal(
      observed_annotations,
      expected_annotations
    )

    opposite_bond <- annotated[
      annotated$protein %in% c("NAT_H_ASPARTATE", "MET_AT_CTERM_BOND"),
      ,
      drop = FALSE
    ]

    expect_identical(
      opposite_bond$specificity,
      rep("semi_Cterm", 2)
    )

    expect_identical(
      opposite_bond$p1_residue,
      c("K", "K")
    )

    expect_identical(
      opposite_bond$p1_prime_residue,
      c("Q", "M")
    )
  }
)
