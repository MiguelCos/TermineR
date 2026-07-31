test_that("diann_adapter reads TSV and Parquet reports without diann", {
  skip_if_not_installed("arrow")

  report <- data.frame(
    Modified.Sequence = c(
      "(UniMod:1)PEPTIDE", "(UniMod:1)PEPTIDE",
      "(UniMod:1)PEPTIDE", "(UniMod:1)PEPTIDE",
      "EPTIDE", "EPTIDE", "EPTIDE",
      "TIDEPE", "TIDEPE", "TIDEPE"
    ),
    Stripped.Sequence = c(
      "PEPTIDE", "PEPTIDE", "PEPTIDE", "PEPTIDE",
      "EPTIDE", "EPTIDE", "EPTIDE",
      "TIDEPE", "TIDEPE", "TIDEPE"
    ),
    Protein.Ids = c(
      "P1", "P1", "P1", "P1", "P2", "P2", "P2", "P3", "P3", "P3"
    ),
    File.Name = c(
      "A.raw", "A.raw", "B.raw", "B.raw", "A.raw", "B.raw", "C.raw",
      "A.raw", "B.raw", "C.raw"
    ),
    Precursor.Normalised = c(10, 5, 20, 3, 8, 15, 18, 4, 7, 9),
    Proteotypic = 1,
    Q.Value = 0.001,
    Protein.Q.Value = 0.001,
    PG.Q.Value = 0.001,
    GG.Q.Value = 0.001,
    check.names = FALSE
  )

  report_dir <- tempfile("terminer-diann-")
  dir.create(report_dir)
  tsv_path <- file.path(report_dir, "report.tsv")
  parquet_path <- file.path(report_dir, "report.parquet")
  data.table::fwrite(report, tsv_path, sep = "\t")
  arrow::write_parquet(report, parquet_path)

  tsv_sum <- suppressWarnings(
    diann_adapter(tsv_path, summarization = "SUM")
  )
  tsv_max <- suppressWarnings(
    diann_adapter(tsv_path, summarization = "MAX")
  )
  parquet_sum <- suppressWarnings(
    diann_adapter(parquet_path, summarization = "SUM")
  )

  expected_columns <- c(
    "nterm_modif_peptide", "nterm_modif", "peptide", "protein"
  )
  expect_true(all(expected_columns %in% names(tsv_sum)))
  expect_true(nrow(tsv_sum) > 0)
  expect_true(nrow(tsv_max) > 0)
  expect_equal(sort(names(tsv_sum)), sort(names(parquet_sum)))
  expect_equal(nrow(tsv_sum), nrow(parquet_sum))
})
