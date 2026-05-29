test_that("UniMod mapping includes Formyl", {
  data("unimod_id_to_name_mapping", package = "TermineR")

  expect_true(any(
    unimod_id_to_name_mapping$id_nr == 122 &
      unimod_id_to_name_mapping$name == "Formyl"
  ))
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
