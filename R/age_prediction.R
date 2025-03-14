
#' Predict Age Score from Gene Expression Matrix
#'
#' Uses pre-trained models to predict age scores based on gene expression.
#'
#' @param mat Matrix. Processed count matrix.
#' @param model_type Character. Model type ('PASTA', 'REG', or 'CT46').
#' @return Numeric vector. Predicted age scores.
#' @export
#' @examples
#' # Assuming 'mat' is your processed expression matrix:
#' # v_age_scores <- predicting_age_score(t(mat))
predicting_age_score <- function(mat, model_type = 'PASTA') {
  requireNamespace('glmnet', quietly = TRUE)
  data(cvfit_REG, envir = environment())
  data(cvfit_PASTA, envir = environment())
  data(beta_PASTA, envir = environment())
  data(cvfit_C46, envir = environment())
  data(beta_C46, envir = environment())
  
  if (model_type == 'PASTA') cur_model <- cvfit_PASTA
  if (model_type == 'REG')   cur_model <- cvfit_REG
  if (model_type == 'CT46')  cur_model <- cvfit_C46
  if (!model_type %in% c('PASTA', 'REG', 'CT46'))
    stop('Specify a valid model; either PASTA, REG, or CT46')
  
  v_age_scores <- stats::predict(cur_model, mat, s = 'lambda.min', 
                                 type = 'link')[, 1] %>% as.numeric
  if (model_type == 'PASTA') v_age_scores <- v_age_scores * beta_PASTA
  if (model_type == 'CT46')  v_age_scores <- v_age_scores * beta_C46
  
  return(v_age_scores)
}

#' Add Age Predictions to Phenotype Data
#'
#' Adds age prediction scores to phenotype data using different models.
#'
#' @param pdata data.table. Phenotype data.
#' @param mat_t Transposed count matrix.
#' @param REG Logical. Include REG model predictions.
#' @param PASTA Logical. Include PASTA model predictions.
#' @param CT46 Logical. Include CT46 model predictions.
#' @return data.table. Updated phenotype data with age predictions.
#' @export
adding_age_preds_to_pdata <- function(pdata, mat_t, REG = TRUE, PASTA = TRUE, CT46 = FALSE) {
  if (REG)   pdata[, REG := predicting_age_score(mat_t, model_type = 'REG')]
  if (PASTA) pdata[, PASTA := predicting_age_score(mat_t, model_type = 'PASTA')]
  if (CT46)  pdata[, CT46  := predicting_age_score(mat_t, model_type = 'CT46')]
  colnames(pdata) %<>% gsub('\\.ch1', '', .)
  colnames(pdata) %<>% gsub('\\:ch1', '', .)
  colnames(pdata) %<>% gsub(' ', '_', .)
  colnames(pdata) %<>% gsub('\\.', '_', .)
  return(pdata)
}

#' Get Phenotype Data with Age Predictions from ExpressionSet
#'
#' Retrieves phenotype data from an ExpressionSet and adds age prediction scores.
#'
#' @param ES ExpressionSet.
#' @param filter_genes Logical. Whether to filter genes using the age model.
#' @param rank_norm Logical. Whether to apply rank normalization.
#' @param REG Logical. Include REG model predictions.
#' @param PASTA Logical. Include PASTA model predictions.
#' @param CT46 Logical. Include CT46 model predictions.
#' @return data.table. Phenotype data with age prediction scores.
#' @export
#' @examples
#' # Assuming you have an ExpressionSet 'ES':
#' # pdata <- getting_pdata_with_age_scores(ES, filter_genes = TRUE, rank_norm = TRUE)
getting_pdata_with_age_scores <- function(ES, filter_genes = FALSE, 
                                          rank_norm = FALSE, REG = TRUE, 
                                          PASTA = TRUE, CT46 = TRUE) {
  mat <- Biobase::exprs(ES)
  if (filter_genes) mat <- filtering_age_model_genes_and_rank_norm(mat)
  if (rank_norm)    mat <- applying_rank_normalization(mat)
  mat_t <- t(mat)
  pdata <- Biobase::pData(ES) %>% data.table::setDT()
  pdata <- adding_age_preds_to_pdata(pdata, mat_t, REG, PASTA, CT46)
  return(pdata)
}
