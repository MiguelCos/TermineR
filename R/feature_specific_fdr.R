#' Apply multiple-testing p-value adjustment on specified features
#'
#' @title feature_fdr_correction
#' @param toptable A data frame with the results of the limma analysis as obtained from the topTable function. It should contain the columns 'P.Value' and 'adj.P.Val' and 'nterm_modif_peptide'.
#' @param interesting_features_table A subset of the peptide data annotation data frame obtained from the `annotate_neo_termini` function. It should contain the list of features that can be defined as 'neo-termini' under the user criteria. Example: All peptides with TMT modification at the N-termini.
#' @param method The method to be used for multiple testing correction. Default is "BH" (Benjamini-Hochberg).
#'
#' @description
#' Apply multiple-testing p-value adjustment on specified features.
#' This approach is based on the idea of independent hypothesis weighting,
#' where the FDR is calculated for a subset of interesting features independently of the tested hypothesis of differential abundance.
#' As described by https://doi.org/10.1038/nmeth.3885
#'
#' @importFrom dplyr right_join mutate filter bind_rows
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#'
#' @export
#' @author Miguel Cosenza-Contreras
feature_fdr_correction <- function (
    toptable,
    interesting_features_table,
    method = "BH") {

  require(dplyr)
  require(magrittr)

  tab_limma_feature_annot <- right_join(
    toptable,
    interesting_features_table
  ) %>%
    mutate(adj.P.Val = p.adjust(p = P.Value, method = method)) %>%
    mutate(fdr_correction = "feature-specific")

  tablimma_subsetout <- filter(toptable, !nterm_modif_peptide %in% interesting_features_table$nterm_modif_peptide) %>%
    mutate(fdr_correction = "global")

  output_limma <- bind_rows(
    tab_limma_feature_annot,
    tablimma_subsetout)

  return(output_limma)
}
