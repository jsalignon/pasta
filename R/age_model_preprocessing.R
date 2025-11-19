
#' Filter Age Model Genes from Count Matrix
#'
#' Subsets the count matrix to include only genes used in the age prediction model.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Filtered count matrix with median imputation.
#' @export
filtering_age_model_genes <- function(mat) {
  data(v_genes_model, envir = environment())
  mat <- mat[match(v_genes_model, rownames(mat)), ]
  median_value <- stats::median(c(mat), na.rm = TRUE)
  mat[is.na(mat)] <- median_value
  rownames(mat) <- v_genes_model
  return(mat)
}

#' Apply Rank Normalization to Matrix
#'
#' Applies rank normalization across each column of the matrix.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Rank-normalized matrix.
#' @export
applying_rank_normalization <- function(mat) {
  mat <- apply(mat, 2, rank, ties.method = 'average')
  return(mat)
}

#' Filter and Rank-Normalize Age Model Genes
#'
#' Combines filtering for age model genes and rank normalization.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Processed count matrix.
#' @export
filtering_age_model_genes_and_rank_norm <- function(mat) {
  mat <- filtering_age_model_genes(mat)
  mat <- applying_rank_normalization(mat)
  return(mat)
}

#' Process GEO Count Matrix for Age Prediction
#'
#' Converts gene IDs and filters & normalizes genes for age prediction.
#'
#' @param mat Matrix. Raw count matrix.
#' @return Matrix. Processed count matrix ready for age prediction.
#' @export
processing_geo_count_mat_for_age_prediction <- function(mat) {
  mat <- converting_entrez_to_ensembl_gene_ids(mat)
  mat <- filtering_age_model_genes_and_rank_norm(mat)
  return(mat)
}

#' Get Processed GEO Count Matrix for Age Prediction
#'
#' Downloads and processes the GEO count matrix for age prediction.
#'
#' @param gse_id Character. GEO series identifier.
#' @return Matrix. Processed count matrix.
#' @export
getting_geo_count_mat_for_age_prediction <- function(gse_id) {
  mat <- getting_geo_count_mat(gse_id)
  mat <- processing_geo_count_mat_for_age_prediction(mat)
  return(mat)
}

#' Convert Mouse Gene Matrix to Human Orthologs
#'
#' Filters a gene expression matrix to retain only mouse genes that have one-to-one 
#' human orthologs, then converts the mouse gene identifiers to their corresponding 
#' human ortholog identifiers. 
#'
#' @param mat Matrix. Gene expression matrix with mouse Ensemble gene identifiers as row
#'   names and samples as columns.
#'
#' @return Matrix. Filtered gene expression matrix with human ortholog identifiers 
#'   as row names, containing only genes with confirmed one-to-one orthology.
#'   
#' @export
keeping_mouse_one2one_orthologs <- function(mat) {
  data(v_human_mouse_one2one, envir = environment())

  mat <- mat[rownames(mat) %in% names(v_human_mouse_one2one), ]
  rownames(mat) <- v_human_mouse_one2one[rownames(mat)]

  return(mat)
}



