
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
#' @param Pasta Logical. Include Pasta model predictions.
#' @param CT46 Logical. Include CT46 model predictions.
#' @return data.table. Phenotype data with age predictions.
#' @export
predicting_age_from_pseudobulks <- function(seu_bulk, 
  REG = TRUE, Pasta = TRUE, CT46 = TRUE){

  requireNamespace('data.table', quietly = TRUE)
  requireNamespace('Seurat', quietly = TRUE)

  mat = Seurat::GetAssayData(seu_bulk, layer = 'data')
  if(inherits(mat, 'dgCMatrix')) { mat %<>% as.matrix }
  mat %<>% filtering_age_model_genes_and_rank_norm
  pdata = seu_bulk[[c('chunk_size', 'type', 'age')]] %>% setDT
  pdata %<>% adding_age_preds_to_pdata(t(mat), REG = REG, Pasta = Pasta, 
    CT46 = CT46)
  return(pdata)
}


#' Aggregate Seurat Data into Pseudobulks and Predict Cellular Age
#'
#' This function aggregates a Seurat object into pseudobulk samples using a specified chunk size and then predicts cellular age scores from the pseudobulk data. It combines the workflow of creating pseudobulk samples via \code{making_pseudobulks_from_seurat} and predicting age using \code{predicting_age_from_pseudobulks}.
#'
#' @param seu A Seurat object containing single-cell transcriptomic data.
#' @param chunk_size Numeric. The number of cells to aggregate into each pseudobulk sample. Default is \code{1000}.
#' @param verbose Logical. If \code{TRUE}, displays progress messages and summaries during processing. Default is \code{TRUE}.
#' @param REG Logical. If \code{TRUE}, include age predictions from the REG model. Default is \code{TRUE}.
#' @param Pasta Logical. If \code{TRUE}, include age predictions from the Pasta model. Default is \code{TRUE}.
#' @param CT46 Logical. If \code{TRUE}, include age predictions from the CT46 model. Default is \code{TRUE}.
#'
#' @return A \code{data.table} containing age prediction scores along with metadata from the pseudobulk samples.
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu' is your Seurat object:
#' pdata <- making_pseudobulks_and_predict_age(seu, chunk_size = 1000)
#' head(pdata)
#' }
#'
#' @export
making_pseudobulks_and_predict_age <- function(seu, chunk_size = 1000, 
  verbose = TRUE, REG = TRUE, Pasta = TRUE, CT46 = TRUE){
  seu_bulk = making_pseudobulks_from_seurat(seu, chunk_size = chunk_size, 
    verbose = verbose)
  pdata = predicting_age_from_pseudobulks(seu_bulk, REG = REG, Pasta = Pasta, 
    CT46 = CT46)
  return(pdata)
}


#' Predict Age Scores Using Multiple Pseudobulk Chunk Sizes
#'
#' This function predicts cellular age scores from a Seurat object by aggregating single-cell data into pseudobulk samples using different chunk sizes. For each specified chunk size, it applies the age prediction workflow and combines the results into a single data.table.
#'
#' @param seu A Seurat object containing single-cell transcriptomic data.
#' @param v_chunk_sizes A numeric vector specifying the different chunk sizes to use when creating pseudobulk samples. Default is \code{c(500, 1000)}.
#' @param REG Logical. If \code{TRUE}, include predictions from the REG model. Default is \code{TRUE}.
#' @param Pasta Logical. If \code{TRUE}, include predictions from the Pasta model. Default is \code{TRUE}.
#' @param CT46 Logical. If \code{TRUE}, include predictions from the CT46 model. Default is \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints progress and summary messages during processing. Default is \code{TRUE}.
#'
#' @return A \code{data.table} containing the combined prediction results from each chunk size. The output includes columns for chunk size, cell type, true age, and predicted age scores for the models used (REG, Pasta, CT46).
#'
#' @examples
#' \dontrun{
#' # Assuming 'seu' is your Seurat object:
#' pdata_big <- predicting_age_multiple_chunks(seu, v_chunk_sizes = c(500, 1000))
#' head(pdata_big)
#' }
#'
#' @export
predicting_age_multiple_chunks <- function(seu, v_chunk_sizes = c(500, 1000), 
  REG = TRUE, Pasta = TRUE, CT46 = TRUE, verbose = TRUE){
  pdata_big = purrr::map(v_chunk_sizes, 
    ~making_pseudobulks_and_predict_age(seu, .x, 
      REG = REG, Pasta = Pasta, CT46 = CT46, verbose = verbose)) %>% 
    do.call(rbind, .)
  return(pdata_big)
}
