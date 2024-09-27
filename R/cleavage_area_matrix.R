#' Calculate positional abundance of amino acid residues from sequences of cleavage areas
#'
#' @title cleavage_area_matrix
#' @param peptides A character vector of amino acid sequences representing cleavage areas of interesting cleavage events. As obtained from the `cleavage_sequence` after annotation with the `annotate_neo_termini` function.
#' @param p_number The number of amino acids to be considered in the cleavage area. Default is 5 for a cleavage area between P5 and P5'.
#'
#' @description
#' Calculate positional abundance of amino acid residues from sequences of cleavage areas.
#' This matrix is a preparation for visualization of the cleavage areas using the pheatmap package.
#'
#' @importFrom dplyr group_by summarize arrange
#' @importFrom purrr map map_dfc
#' @importFrom stringr str_split
#' @importFrom tidyr pivot_longer pivot_wider complete
#' @importFrom tibble as_tibble column_to_rownames tibble
#' @importFrom magrittr %>%
#'
#' @export
#' @author Miguel Cosenza-Contreras
cleavage_area_matrix <- function(peptides,
                                 p_number = 5) {
  
  require(dplyr)
  require(purrr)
  require(stringr)
  require(tidyr)
  require(tibble)
  require(magrittr)
  
  # Split each peptide into individual amino acids using map and str_split
  split_peptides <- purrr::map(peptides, str_split, "")
  
  # Convert the list of split peptides into a data frame using map_dfc
  peptide_matrix <- map_dfc(split_peptides, as.data.frame) %>%
    # Transpose the data frame so that each row corresponds to a peptide
    t() %>%
    suppressMessages()
  
  # Create character vector that goes from PX to P1 and then P1' to PX'
  positions <- c(paste0("P", p_number:1),
                 paste0("P", 1:p_number, "'"))
  
  # Define the order of amino acids for sorting
  amino_acid_order <- c("A", "C", "D", "E", "F", "G", "H", 
                        "I", "K", "L", "M", "N", "P", "Q", 
                        "R", "S", "T", "V", "W", "Y")
  
  # Use the positions object to define the column names of the peptide sequence matrix
  colnames(peptide_matrix) <- positions
  rownames(peptide_matrix) <- NULL
  
  # Count how many times a particular amino acid appears in a certain position
  amino_acid_count <- peptide_matrix %>%
    as_tibble() %>%
    pivot_longer(cols = everything(),
                 names_to = "position",
                 values_to = "AA") %>%
    group_by(position, AA) %>%
    summarize(n = n(), .groups = "drop") %>%
    complete(position, AA = amino_acid_order, fill = list(n = 0)) %>%
    pivot_wider(id_cols = AA, names_from = position, values_from = n)
  
  # Arrange the rows based on the order of the amino acids in the vector above
  amino_acid_count <- amino_acid_count %>%
    arrange(factor(AA, levels = amino_acid_order)) %>%
    column_to_rownames("AA")
  
  # Select only columns between P5 and P5' from the count matrix
  amino_acid_count <- amino_acid_count[,
                                       c(paste0("P", p_number:1),
                                         paste0("P", 1:p_number, "'"))] %>%
    as.matrix()
  
  # Return the final matrices as a list
  return(list(peptide_matrix = peptide_matrix,
              amino_acid_count = amino_acid_count))
}

