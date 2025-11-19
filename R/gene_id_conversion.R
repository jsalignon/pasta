
#' Convert Entrez to Ensembl Gene IDs (using GEO names)
#'
#' Replaces row names in the matrix with Ensembl gene IDs from a GEO-specific mapping.
#'
#' @param mat Matrix. Count matrix with original gene identifiers.
#' @return Matrix. Count matrix with Ensembl gene IDs as row names.
#' @export
converting_entrez_to_ensembl_gene_ids <- function(mat) {
  data(v_new_names_geo_mat, envir = environment())
  rownames(mat) <- v_new_names_geo_mat
  mat <- mat[grep('ENSG', v_new_names_geo_mat, value = TRUE), ]
  return(mat)
}

#' Convert Gene IDs to Ensembl Gene IDs using biomaRt (Placeholder)
#'
#' This function is a placeholder for gene conversion using biomaRt.
#'
#' @param mat Matrix. Count matrix.
#' @return Matrix. Count matrix with Ensembl gene IDs.
#' @export
converting_gene_ids_to_ensembl_gene_ids_with_biomaRt <- function(mat) {
  rownames(mat) <- v_new_names_geo_mat
  mat <- mat[grep('ENSG', v_new_names_geo_mat, value = TRUE), ]
  return(mat)
}

#' Get Ensembl Human Mart
#'
#' Returns a biomaRt Mart object for human gene annotation.
#'
#' @return Mart object.
#' @export
get_mart_ensembl_human <- function() {
  biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
}

#' Check biomaRt Attributes
#'
#' Lists attributes from a biomaRt Mart filtered by a pattern.
#'
#' @param cur_mart Mart object.
#' @param pattern Character. Regular expression to filter attributes.
#' @return Printed data frame of matching attributes.
#' @export
check_biomart_attributes <- function(cur_mart, pattern) {
  biomaRt::listAttributes(cur_mart) %>% 
    .[grep(pattern, .[,1], ignore.case = TRUE), ] %>% 
    print
}

#' Convert Gene IDs to Ensembl Gene IDs using biomaRt
#'
#' Uses biomaRt to map a platform-specific gene identifier to Ensembl gene IDs.
#'
#' @param mat Matrix. Count matrix with row names as gene IDs.
#' @param cur_mart Mart object.
#' @param platform Character. Platform identifier (e.g. 'illumina_humanht_12_v4').
#' @return data.table. Mapping between platform IDs and Ensembl gene IDs.
#' @export
converting_gene_ids_to_ensembl_gene_ids <- function(mat, cur_mart,
                                                   platform = 'illumina_humanht_12_v4') {
  requireNamespace('biomaRt', quietly = TRUE)
  v_ids <- rownames(mat)
  dt_ids <- biomaRt::getBM(
    attributes = c(platform, 'ensembl_gene_id'),
    filters    = platform,
    values     = v_ids,
    mart       = cur_mart
  ) %>% data.table::setDT()
  return(dt_ids)
}

#' Rename Matrix Rows after Gene ID Conversion
#'
#' Renames the rows of a matrix based on a conversion table.
#'
#' @param mat Matrix. Count matrix.
#' @param dt_ids data.table. Conversion table with original and Ensembl gene IDs.
#' @param species character. Either mouse or human.
#' @return Matrix. Count matrix with new row names.
#' @export
renaming_rows_in_mat_after_gene_id_conversion <- function(mat, dt_ids, 
  species = c('human','mouse')) {
  species <- match.arg(species)
  v_names <- stats::setNames(dt_ids$ensembl_gene_id, dt_ids[[1]])
  v_names[duplicated(v_names)] <- NA
  v_new_names <- v_names[rownames(mat)]
  v_new_names[is.na(v_new_names)] <- paste0('NA_', seq_along(which(is.na(v_new_names))))
  names(v_new_names) <- NULL
  mat1 <- mat
  rownames(mat1) <- v_new_names
  # mat1 <- mat1[grep('ENSG', v_new_names, value = TRUE), ]

  if(species == 'human'){
    mat1 <- mat1[grep('ENSG', v_new_names, value = TRUE), ]
  } else if (species == 'mouse'){
    mat1 <- mat1[grep('ENSMUSG', v_new_names, value = TRUE), ]
  }

  return(mat1)
}
