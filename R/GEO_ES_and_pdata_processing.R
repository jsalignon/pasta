
#' Retrieve GEO Phenotype Data
#'
#' Retrieves and combines phenotype data from GEO.
#'
#' @param gse_id Character. GEO series identifier.
#' @param mat Matrix. Count matrix (not used here but maintained for consistency).
#' @return data.frame. Phenotype data.
#' @export
getting_GEO_pdata <- function(gse_id, mat) {
  lES <- GEOquery::getGEO(gse_id, getGPL = FALSE)
  if (length(lES) > 1) {
    warning('Multiple data sources found. Combining them may fail or give incorrect results. Please double-check.')
    l_pdata <- purrr::map(lES, ~Biobase::pData(.x))
    raw_pdata <- do.call(rbind, l_pdata)
    rownames(raw_pdata) <- raw_pdata$geo_accession
  } else {
    ES <- lES[[1]]
    raw_pdata <- data.frame(Biobase::pData(ES))
  }
  return(raw_pdata)
}

#' Create ExpressionSet from GEO Data
#'
#' Constructs an ExpressionSet object from a count matrix and phenotype data.
#'
#' @param mat Matrix. Count matrix.
#' @param pdata data.frame. Phenotype data.
#' @return ExpressionSet. ExpressionSet object.
#' @export
getting_GEO_ES <- function(mat, pdata) {
  samples_mat <- colnames(mat)
  samples_pda <- rownames(pdata)
  if (!(length(samples_mat) == length(samples_pda)))
    stop('Sample length differs between the count and annotation matrices')
  if (!all.equal(samples_mat, samples_pda))
    stop('Sample IDs differ between the count and annotation matrices')
  ES <- Biobase::ExpressionSet(mat, phenoData = Biobase::AnnotatedDataFrame(pdata))
  return(ES)
}

#' Get GEO ExpressionSet for Age Model
#'
#' Downloads the GEO count matrix and phenotype data and returns an ExpressionSet.
#'
#' @param gse_id Character. GEO series identifier.
#' @return ExpressionSet.
#' @export
getting_GEO_ES_for_age_model <- function(gse_id) {
  mat <- getting_geo_count_mat_for_age_prediction(gse_id)
  pdata <- getting_GEO_pdata(gse_id, mat)
  sel_samples <- pdata$geo_accession %>% .[. %in% colnames(mat)]
  mat1 <- mat[, sel_samples]
  pdata1 <- pdata %>% .[match(sel_samples, .$geo_accession), ]
  len_pdata <- nrow(pdata)
  len_mat <- length(sel_samples)
  if (len_mat != len_pdata)
    print(paste(len_pdata - len_mat, 'samples missing in the matrix'))
  ES <- getting_GEO_ES(mat1, pdata1)
  return(ES)
}
