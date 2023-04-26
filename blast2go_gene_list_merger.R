#' Merge Blast2GO Gene List
#'
#' This function takes the output of a Blast2GO analysis and a gene list file,
#' and merges them based on the gene ID. The resulting file will include the
#' gene description, sequence length, E-value, hit name, hit description, and
#' Gene Ontology (GO) terms for each gene.
#'
#' @param blast2go_file Path to the Blast2GO output file.
#' @param blast2go_columns Column names of the Blast2GO file as a character vector.
#' @param gene_list_file Path to the gene list file.
#' @param merged_file Path to the file to save the merged data to.
#' @param remove_endings Logical indicating whether to remove non-alphanumeric endings from gene IDs in the Blast2GO file. Default is TRUE.
#'
#' @return A TSV file containing the merged data.
#'
#' @importFrom utils write.table
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unnest
#' @importFrom stringr str_remove
#' @importFrom magrittr %>%
#' @importFrom readr read_delim
#' @importFrom purrr possibly
#'
#' @author Carson Stacy <carsonstacy2@gmail.com>
#' @export
#'
#' @examples
#' merge_blast2go_gene_list("blast2go_output.txt", c("Gene_ID", "Sequence_Length", "Sequence_Name", "E_Value", "Hit_ID", "Hit_Name", "Hit_Description", "GO"), "gene_list.txt", "merged_data.txt")
#'
#' @seealso https://www.blast2go.com/
#' @seealso https://en.wikipedia.org/wiki/Gene_ontology
merge_blast2go_gene_list <- function(blast2go_file, blast2go_columns, gene_list_file, merged_file, remove_endings = TRUE) {
  # Check if blast2go_columns is a character vector
  if (!is.character(blast2go_columns)) {
    stop("blast2go_columns must be a character vector.")
  }
  
  # Read in the Blast2GO output file without column names
  blast2go <- tryCatch({
    read.delim(blast2go_file, header = FALSE, sep = "\t")
  }, error = function(e) {
    stop("Error reading Blast2GO file. Check if the file path is correct and the file is in the correct format.")
  })
  
  # Add column names to the Blast2GO data frame
  colnames(blast2go) <- blast2go_columns
  
  # Check if remove_endings is a logical value
  if (!is.logical(remove_endings)) {
    stop("remove_endings must be a logical value.")
  }
  
    # Remove any non-alphanumeric endings from the gene IDs in the Blast2GO output file
  if (remove_endings) {
    blast2go$Gene_ID <- tryCatch(stringr::str_remove(blast2go$Gene_ID, "\\W+$"), error = function(e) {
      stop("Error removing non-alphanumeric endings from Blast2GO gene IDs: ", e$message)
    })
  }
  
  # Group the Blast2GO data frame by gene ID and concatenate the GO terms using "|"
  blast2go_grouped <- tryCatch(aggregate(GO, by = list(Gene_ID = blast2go$Gene_ID), FUN = paste, collapse = "|"), error = function(e) {
    stop("Error grouping Blast2GO data frame: ", e$message)
  })
  
  # Join the gene description column to the grouped data frame
  blast2go_grouped <- tryCatch(merge(blast2go_grouped, blast2go[c("Gene_ID", "Hit_Description")], by = "Gene_ID"), error = function(e) {
    stop("Error joining gene description column to Blast2GO data frame: ", e$message)
  })
  
  # Rename the columns in the merged data frame
  tryCatch(colnames(blast2go_grouped) <- c("Gene_ID", "GO_terms", "Gene_Description"), error = function(e) {
    stop("Error renaming columns in merged data frame: ", e$message)
  })
  
  # Read in the gene list file
  gene_list <- tryCatch(read.delim(gene_list_file, header = TRUE, sep = "\t"), error = function(e) {
    stop("Error reading gene list file: ", e$message)
  })
  
  # Merge the grouped data frame with the gene list
  merged_data <- tryCatch(merge(gene_list, blast2go_grouped, by = "Gene_ID", all.x = TRUE), error = function(e) {
    stop("Error merging data frames: ", e$message)
  })
  
  # Write the merged data frame to a new TSV file
  tryCatch(write.table(merged_data, merged_file, sep = "\t", quote = FALSE, row.names = FALSE), error = function(e) {
    stop("Error writing merged data frame to file: ", e$message)
  })
  
  # Print a message to indicate completion of the function
  message("Merged data saved to", merged_file)
}
