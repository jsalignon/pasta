
#' Check if URL is valid
#'
#' Attempts to open a connection to a URL to check its validity.
#'
#' @param url_in Character. The URL to be checked.
#' @param t Numeric. Timeout in seconds.
#' @return Logical. TRUE if the URL is valid, FALSE otherwise.
#' @keywords internal
checking_if_url_is_valid <- function(url_in, t = 2) {
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = 'rt', timeout = t),
                                  silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))
  ifelse(is.null(check), TRUE, FALSE)
}

#' Construct GEO Download URL
#'
#' Constructs the download URL for GEO count matrix files.
#'
#' @param gse_id Character. GEO series identifier.
#' @return Character. The URL to download the count matrix.
#' @keywords internal
getting_geo_url <- function(gse_id) {
  base_url <- 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts'
  geo_url <- paste0(base_url, '&acc=', gse_id, '&file=', gse_id, 
                    '_raw_counts_GRCh38.p13_NCBI.tsv.gz')
  return(geo_url)
}

#' Validate GEO URL and Return It
#'
#' Checks if the GEO URL is valid and either returns the URL or stops execution.
#'
#' @param gse_id Character. GEO series identifier.
#' @param return_boolean Logical. If TRUE, only return the validation status.
#' @return Character or Logical. The URL if valid, or validation status if requested.
#' @keywords internal
getting_geo_valid_url <- function(gse_id, return_boolean = FALSE) {
  geo_url <- getting_geo_url(gse_id)
  valid_url <- checking_if_url_is_valid(geo_url)
  if (return_boolean) return(valid_url)
  if (!valid_url) stop('Count matrix not available for this GSE id.')
  return(geo_url)
}

#' Download GEO Count Matrix
#'
#' Downloads the GEO count matrix and returns it as a matrix.
#'
#' @param gse_id Character. GEO series identifier.
#' @return Matrix. Count matrix with integer values.
#' @export
getting_geo_count_mat <- function(gse_id) {
  geo_url <- getting_geo_valid_url(gse_id)
  mat <- as.matrix(data.table::fread(geo_url, header = TRUE, colClasses = 'integer'),
                   rownames = 1)
  return(mat)
}
