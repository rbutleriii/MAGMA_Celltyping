#' Standardise the column headers in the Summary Statistics files
#'
#' Use a reference data table of common column header names (stored in sumstatsColHeaders.rda) convert them to a standard set, i.e. chromosome --> CHR
#' 
#' This function does not check that all the required column headers are present.
#' 
#' The amended header is returned as a character, and is not written directly back into the file.
#'
#' @param first_line Column header (character)
#'
#' @return The amended column header (the column headers will not be written directly into the summary statistics file)
#'
#' @examples
#' col_headers = standardise.sumstats.column.headers.crossplatform(first_line)
#'
#' @export

standardise.sumstats.column.headers.crossplatform <- function (first_line) 
{
  column_headers = strsplit(first_line, "\t")[[1]]
  print("First line of summary statistics file: ")
  print(first_line)
  print(column_headers)
  data(sumstatsColHeaders) # Loads the pre-defined "common column names" for GWAS data
  
  column_headers = toupper(column_headers)
  print(column_headers)
  
  for (headerI in 1:dim(sumstatsColHeaders)[1]) {
    un = sumstatsColHeaders[headerI, 1]
    cr = sumstatsColHeaders[headerI, 2]
    if (un %in% column_headers & (!cr %in% column_headers)) {
      column_headers = gsub(sprintf("^%s$", un), cr, column_headers)
    }
  }
  new_first_line = paste(column_headers, collapse = "\t")
  return(new_first_line)
}
