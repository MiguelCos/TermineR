test_that("UniMod mapping includes Formyl", {
  data("unimod_id_to_name_mapping", package = "TermineR")

  expect_true(all(
    c(
      "id",
      "id_nr",
      "name",
      "monoisotopic_mass",
      "average_mass",
      "composition"
    ) %in% colnames(unimod_id_to_name_mapping)
  ))

  formyl_mapping <- unimod_id_to_name_mapping[
    unimod_id_to_name_mapping$id_nr == 122,
  ]

  expect_true(any(
    unimod_id_to_name_mapping$id_nr == 122 &
      unimod_id_to_name_mapping$name == "Formyl"
  ))
  expect_equal(formyl_mapping$monoisotopic_mass[[1]], 27.994915, tolerance = 1e-6)
  expect_equal(formyl_mapping$average_mass[[1]], 28.0101, tolerance = 1e-6)
  expect_equal(formyl_mapping$composition[[1]], "C O")
})

test_that("TargetP package data have expected columns", {
  targetp_objects <- c(
    "human_targetp_processing",
    "mouse_targetp_processing",
    "arabidopsis_targetp_processing",
    "rat_targetp_processing",
    "yeast_targetp_processing",
    "medicago_trucantula_targetp_processing",
    "rhizobium_melitoli_targetp_processing",
    "pig_targetp_processing",
    "human_iso_targetp_processing",
    "ecoli_targetp_processing",
    "c_elegans_targetp_processing",
    "synechocystis_targetp_processing"
  )

  for(targetp_object in targetp_objects){

    targetp_env <- new.env(parent = emptyenv())
    data(list = targetp_object, package = "TermineR", envir = targetp_env)
    targetp_processing <- get(targetp_object, envir = targetp_env)

    expect_true(all(
      c("protein", "targetp_category", "targetp_p1_position") %in%
        colnames(targetp_processing)
    ))
    expect_true(nrow(targetp_processing) > 0)
    expect_true(all(!is.na(targetp_processing$targetp_p1_position)))

  }
})

test_that("CaspSites package data have expected columns", {
  caspsites_env <- new.env(parent = emptyenv())
  data("human_caspsites_processing", package = "TermineR", envir = caspsites_env)
  human_caspsites_processing <- get("human_caspsites_processing", envir = caspsites_env)

  expect_true(all(
    c(
      "protein",
      "caspsites_p1_position",
      "caspsites_p1_prime_position",
      "caspsites_cleavage_sequence",
      "caspsites_cleavage_site",
      "caspsites_datasets",
      "caspsites_evidence_count"
    ) %in% colnames(human_caspsites_processing)
  ))

  expect_true(nrow(human_caspsites_processing) > 0)
  expect_true(all(!is.na(human_caspsites_processing$protein)))
  expect_true(all(!is.na(human_caspsites_processing$caspsites_p1_prime_position)))
  expect_true(all(human_caspsites_processing$caspsites_p1_prime_position > 0))
  caspsites_sequences <- unlist(strsplit(
    human_caspsites_processing$caspsites_cleavage_sequence,
    "|",
    fixed = TRUE
  ))
  expect_true(all(nchar(caspsites_sequences) == 8))
})
test_that("annotate_neo_termini adds human CaspSites exact-site annotation", {
  peptide <- "GSASEVPSELSERPK"
  sequence <- paste0(strrep("A", 40), "DQKD", peptide, strrep("A", 20))
  fasta_file <- tempfile(fileext = ".fasta")

  writeLines(
    c(
      ">sp|A0A096LP01|SIM26_HUMAN",
      sequence
    ),
    fasta_file
  )

  peptides_df <- data.frame(
    nterm_modif_peptide = paste0("n_", peptide),
    nterm_modif = "n",
    peptide = peptide,
    protein = "A0A096LP01",
    sample1 = 1,
    check.names = FALSE
  )

  annotated <- annotate_neo_termini(
    peptides_df = peptides_df,
    fasta_location = fasta_file,
    sense = "C",
    specificity = "K|R",
    organism = "human",
    pssm_prediction = FALSE
  )

  expect_true(any(annotated$caspsites_matches_p1_prime))
  expect_true(any(annotated$caspsites_p1_prime_position == 45))
  expect_true(any(grepl("Caspase-3", annotated$caspsites_datasets, fixed = TRUE)))
})
test_that("unsupported organisms fail early", {
  expect_error(
    annotate_neo_termini(
      peptides_df = data.frame(),
      fasta_location = "missing.fasta",
      sense = "C",
      specificity = "K|R",
      organism = "dog"
    ),
    "Supported organisms"
  )
})

test_that("FragPipe label-free adapter maps N-terminal Formyl", {
  parent_dir <- tempfile("fragpipe_lf")
  dir.create(parent_dir)

  psm_file <- file.path(parent_dir, "psm.tsv")
  annotation_file <- tempfile("annotation", fileext = ".txt")

  readr::write_tsv(
    data.frame(
      `Spectrum File` = "interact-run1.pep.xml",
      Peptide = "PEPTIDE",
      `Modified Peptide` = "n[27.994915]PEPTIDE",
      Probability = 0.99,
      Intensity = 1000,
      `Assigned Modifications` = "N-term(27.994915)",
      `Is Unique` = TRUE,
      `Protein ID` = "P12345",
      Gene = "GENE1",
      check.names = FALSE
    ),
    psm_file
  )

  readr::write_tsv(
    data.frame(
      run = "run1",
      sample = "sample1"
    ),
    annotation_file
  )

  result <- fragpipe_lf_adapter(
    parent_dir = parent_dir,
    annotation_file_path = annotation_file
  )

  expect_true("Formyl_PEPTIDE" %in% result$nterm_modif_peptide)
  expect_equal(result$nterm_modif[[1]], "Formyl")
})
