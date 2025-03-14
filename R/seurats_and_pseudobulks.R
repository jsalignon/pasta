
#' Filter Cell Types in a Seurat Object
#'
#' Filters a Seurat object to include only cell types with a minimum number of cells.
#'
#' @param seu Seurat object.
#' @param n_cell_min Numeric. Minimum number of cells required per cell type.
#' @param dry_run Logical. If TRUE, do not return filtered object.
#' @param verbose Logical. If TRUE, prints cell type summaries.
#' @return Seurat object (if dry_run is FALSE).
#' @export
filter_cell_types_in_seu_object <- function(seu, n_cell_min = 500, 
                                              dry_run = FALSE, verbose = FALSE) {
  v_type_to_keep <- seu$type %>% table %>% .[. >= n_cell_min] %>% names
  seu1 <- seu[, seu$type %in% v_type_to_keep]
  if (verbose) print(table(seu1[[c('type', 'age')]]))
  if (!dry_run) return(seu1)
}

#' Create Pseudobulk Samples from Seurat Object
#'
#' Aggregates single-cell data into pseudobulk samples based on cell type and age.
#'
#' @param seu Seurat object.
#' @param chunk_size Numeric. Number of cells per pseudobulk sample.
#' @param verbose Logical. If TRUE, prints summary tables.
#' @return Seurat object with pseudobulk aggregated data.
#' @export
making_pseudobulks_from_seurat <- function(seu, chunk_size = 1000, verbose = TRUE) {
  requireNamespace('data.table', quietly = TRUE)
  requireNamespace('Seurat', quietly = TRUE)
  
  if (chunk_size == 1) {
    seu_bulk <- seu
  } else {
    pdata1 <- seu[[c('type', 'age')]] %>% data.table::setDT()
    pdata1[, id := 1:.N, by = c('type', 'age')]
    pdata1[, chunk := round(id / chunk_size) + 1]
    pdata1[, type_1 := paste0(type, '_', chunk)]
    seu$type_1 <- pdata1$type_1
    
    n_types <- length(unique(seu$type))
    if (n_types == 1) {
      v_grouping_vars <- c('type_1', 'age')
    } else {
      v_grouping_vars <- c('type_1', 'age', 'type')
    }
    
    seu_bulk <- Seurat::AggregateExpression(seu, 
                                            group.by = v_grouping_vars, 
                                            return.seurat = TRUE, verbose = FALSE)
    
    if (n_types == 1) seu_bulk$type <- unique(seu$type)
  }
  
  seu_bulk$age <- as.numeric(seu_bulk$age)
  seu_bulk$chunk_size <- chunk_size
  
  if (verbose) seu_bulk[[c('age', 'type')]] %>% table %>% t %>% print
  return(seu_bulk)
}

#' Predict Age from Pseudobulk Data
#'
#' Uses processed pseudobulk data to predict age scores.
#'
#' @param seu_bulk Seurat object with pseudobulk data.
#' @param REG Logical. Include REG model predictions.
#' @param PASTA Logical. Include PASTA model predictions.
#' @param CT46 Logical. Include CT46 model predictions.
#' @return data.table. Phenotype data with age predictions.
#' @export
predicting_age_from_pseudobulks <- function(seu_bulk, 
  REG = TRUE, PASTA = TRUE, CT46 = TRUE){

  requireNamespace('data.table', quietly = TRUE)
  requireNamespace('Seurat', quietly = TRUE)

  mat = Seurat::GetAssayData(seu_bulk, layer = 'data')
  if(inherits(mat, 'dgCMatrix')) { mat %<>% as.matrix }
  mat %<>% filtering_age_model_genes_and_rank_norm
  pdata = seu_bulk[[c('chunk_size', 'type', 'age')]] %>% setDT
  pdata %<>% adding_age_preds_to_pdata(t(mat), REG = REG, PASTA = PASTA, 
    CT46 = CT46)
  return(pdata)
}
