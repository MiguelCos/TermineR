test_that("protein adjustment returns ratios and preserves signed log2 values", {
  peptides <- data.frame(
    nterm_modif_peptide = c("tryptic_a", "tryptic_b", "target", paste0("background_", 1:7), "low"),
    protein = c(rep("P_target", 3), paste0("P_background_", 1:7), "P_low"),
    sample_a.intensity = c(50, 50, 25, rep(64, 7), 0.25),
    sample_b.intensity = c(100, 100, 50, rep(64, 7), 0.25),
    sample_c.intensity = c(50, 50, 50, rep(64, 7), 0.25)
  )
  annot <- data.frame(sample = c("sample_a", "sample_b", "sample_c"))
  peptide_annot <- data.frame(
    nterm_modif_peptide = peptides$nterm_modif_peptide,
    specificity = c("specific", "specific", "semi_Nterm", rep("specific", 8))
  )

  result <- suppressMessages(peptide2protein_normalization(peptides, annot, peptide_annot))
  ratios <- result$protein_normalized_pepts_abundance
  target <- ratios[ratios$nterm_modif_peptide == "target", -1]
  expect_equal(as.numeric(target[1, ]), c(0.25, 0.25, 0.5))
  expect_equal(target[[2]] / target[[1]], 1)
  expect_equal(target[[3]] / target[[1]], 2)
  expect_identical(names(ratios), c("nterm_modif_peptide", paste0("fraction_int_peptide2prot_", annot$sample)))

  # Stable background peptides keep sample and global medians equal.
  scaled <- result$protein_normalized_pepts_scaled
  expect_equal(as.numeric(scaled[scaled$nterm_modif_peptide == "target", -1]), c(-2, -2, -1))
  proteins <- result$summarized_protein_abundance
  expect_equal(as.numeric(proteins[proteins$protein == "P_target", -1]), c(100, 200, 100))
  expect_identical(names(proteins), c("protein", paste0(annot$sample, ".intensity_prot")))
  scaled_proteins <- result$summarized_protein_abundance_scaled
  expect_equal(as.numeric(scaled_proteins[scaled_proteins$protein == "P_low", -1]), rep(-2, 3))
  expect_identical(names(result), c(
    "protein_normalized_pepts_scaled", "protein_normalized_pepts_abundance",
    "summarized_protein_abundance", "summarized_protein_abundance_scaled",
    "summarize_by_specificity"
  ))
  expect_true(result$summarize_by_specificity)

  # The existing option to summarize all peptides still controls the denominator.
  all_peptides <- suppressMessages(peptide2protein_normalization(
    peptides, annot, peptide_annot, summarize_by_specificity = FALSE
  ))
  all_ratios <- all_peptides$protein_normalized_pepts_abundance
  expect_equal(as.numeric(all_ratios[all_ratios$nterm_modif_peptide == "target", -1]),
               c(25 / 125, 50 / 250, 50 / 150))
  expect_false(all_peptides$summarize_by_specificity)
})

test_that("invalid protein denominators yield NA and retain downstream row handling", {
  peptides <- data.frame(
    nterm_modif_peptide = c(paste0("tryptic_", 1:6), paste0("target_", 1:7), paste0("background_", 1:3)),
    protein = c(paste0("P", 1:6), paste0("P", 1:7), paste0("BG", 1:3)),
    sample_a = c(0, -1, Inf, -Inf, NA_real_, NaN, rep(8, 7), rep(64, 3))
  )
  peptides$sample_b <- peptides$sample_a
  annot <- data.frame(sample = c("sample_a", "sample_b"))
  peptide_annot <- data.frame(
    nterm_modif_peptide = peptides$nterm_modif_peptide,
    specificity = c(rep("specific", 6), rep("semi_Nterm", 7), rep("specific", 3))
  )

  result <- suppressMessages(peptide2protein_normalization(peptides, annot, peptide_annot))
  ratios <- result$protein_normalized_pepts_abundance
  expect_true(all(is.na(ratios[ratios$nterm_modif_peptide %in% paste0("target_", 1:7), -1])))
  expect_false(any(is.infinite(as.matrix(ratios[, -1]))))
  expect_identical(result$protein_normalized_pepts_scaled$nterm_modif_peptide, paste0("background_", 1:3))
  expect_equal(as.matrix(result$protein_normalized_pepts_scaled[, -1]), matrix(
    0, nrow = 3, ncol = 2,
    dimnames = list(NULL, paste0("fraction_int_peptide2prot_", annot$sample))
  ), ignore_attr = TRUE)

  # Protein summaries stay raw; only the matrix used for log2 is masked.
  proteins <- result$summarized_protein_abundance
  expect_equal(proteins$sample_a_prot[match(c("P1", "P2", "P5", "P6"), proteins$protein)], c(0, -1, 0, 0))
  expect_equal(proteins$sample_a_prot[match(c("P3", "P4"), proteins$protein)], c(Inf, -Inf))
  expect_false("P7" %in% proteins$protein)
  scaled_proteins <- result$summarized_protein_abundance_scaled
  expect_true(all(is.na(scaled_proteins[scaled_proteins$protein %in% paste0("P", 1:6), -1])))
  expect_false(any(is.infinite(as.matrix(scaled_proteins[, -1]))))
  expect_equal(scaled_proteins$sample_a_prot[match(paste0("BG", 1:3), scaled_proteins$protein)], rep(6, 3))
})
