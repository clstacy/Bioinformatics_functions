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
#'
#' @author Carson Stacy <carsonstacy2@gmail.com>
#' @export
#'
#' @examples
#' merge_blast2go_gene_list("blast2go_output.txt", c("Gene_ID", "Sequence_Length", "Sequence_Name", "E_Value", "Hit_ID", "Hit_Name", "Hit_Description", "GO"), "gene_list.txt", "merged_data.txt")
#'
#' @seealso https://www.blast2go.com/
#' @seealso https://en.wikipedia.org/wiki/Gene_ontology
merge_blast2go_gene_list <- function(blast2go_file, blast2go_columns, gene_list_file, merged_file) {
  # Read in the Blast2GO output file without column names
  blast2go <- read.delim(blast2go_file, header = FALSE, sep = "\t")
  
  # Add column names to the Blast2GO data frame
  colnames(blast2go) <- blast2go_columns
  
  # Remove any non-alphanumeric endings from the gene IDs in the Blast2GO output file
  blast2go$Gene_ID <- stringr::str_remove(blast2go$Gene_ID, "\\W+$")
  
  # Group the Blast2GO data frame by gene ID and concatenate the GO terms using "|"
  blast2go_grouped <- aggregate(GO, by = list(Gene_ID = blast2go$Gene_ID), FUN = paste, collapse = "|")
  
  # Join the gene description column to the grouped data frame
  blast2go_grouped <- merge(blast2go_grouped, blast2go[c("Gene_ID", "Hit_Description")], by = "Gene_ID")
  
  # Rename the columns in the merged data frame
  colnames(blast2go_grouped) <- c("Gene_ID", "GO_terms", "Gene_Description")
  
  # Read in the gene list file
  gene_list <- read.delim(gene_list_file, header = TRUE, sep = "\t")
  
  # Merge the grouped data frame with the gene list
  merged_data <- merge(gene_list, blast2go_grouped, by = "Gene_ID", all.x = TRUE)
  
  # Write the merged data frame to a new TSV file
  write.table(merged_data, merged_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Print a message to indicate completion of the function
  cat("Merged data saved to", merged_file, "\n")
}
