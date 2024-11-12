
utils::globalVariables(c('v_genes_model', 'v_new_names_geo_mat', 'cvfit_REG', 
	'cvfit_PASTA', 'beta_PASTA', 'cvfit_C46', 'beta_C46', 'ES_GSE103938',
	'seu_orozco_2020_retina_horizontal_cells'))


# magrittr and data.table variables
utils::globalVariables(c('.', 'chunk_size', 'id', 'chunk', 'type_1', 'type'))


#' Title of the function
#'
#' Description of the function.
#'
#' @export
#' @keywords internal
checking_if_url_is_valid <- function(url_in,t=2){
  con = url(url_in)
  check = suppressWarnings(try(open.connection(con, open = 'rt', timeout = t),
  	silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))
  ifelse(is.null(check), TRUE, FALSE)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @keywords internal
getting_geo_url <- function(gse_id){
	urld    = 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts'
	geo_url = paste0(urld, '&acc=', gse_id, '&file=', gse_id, 
		'_raw_counts_GRCh38.p13_NCBI.tsv.gz')
	return(geo_url)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @keywords internal
getting_geo_valid_url <- function(gse_id){
	geo_url = getting_geo_url(gse_id)
	valid_url = checking_if_url_is_valid(geo_url)
  if(!valid_url) stop('Count matrix not available for this GSE id.')
	return(geo_url)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
getting_geo_count_mat <- function(gse_id){
	urld <- 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts'
	path <- paste0(urld, '&acc=', gse_id, '&file=', gse_id, '_raw_counts_GRCh38.p13_NCBI.tsv.gz');
	geo_url = getting_geo_valid_url(gse_id)
	mat <- as.matrix(data.table::fread(geo_url, header = TRUE, 
		colClasses = 'integer'), rownames = 1)
	return(mat)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
converting_entrez_to_ensembl_gene_ids <- function(mat){
	# Loading the PASTA model
  data(v_new_names_geo_mat, envir = environment())

	rownames(mat) = v_new_names_geo_mat
	mat = mat[grep('ENSG', v_new_names_geo_mat, value = TRUE),]
	return(mat)

}


#' Title of the function
#'
#' Description of the function.
#'
#' @export
filtering_age_model_genes <- function(mat){
	data(v_genes_model, envir = environment())
	mat = mat[match(v_genes_model, rownames(mat)), ]
	median_value = stats::median(c(mat), na.rm = TRUE)
	mat[is.na(mat)] = median_value
	rownames(mat) = v_genes_model
	return(mat)
}


#' Title of the function
#'
#' Description of the function.
#'
#' @export
applying_rank_normalization <- function(mat){
	mat = apply(mat, 2, rank, ties.method = 'average')
	return(mat)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
filtering_age_model_genes_and_rank_norm <- function(mat){
	mat %<>% filtering_age_model_genes
	mat %<>% applying_rank_normalization
	return(mat)
}


#' Title of the function
#'
#' Description of the function.
#'
#' @export
processing_geo_count_mat_for_age_prediction <- function(mat){
	mat %<>% converting_entrez_to_ensembl_gene_ids
	mat %<>% filtering_age_model_genes_and_rank_norm
	return(mat)
}


#' Title of the function
#'
#' Description of the function.
#'
#' @export
getting_geo_count_mat_for_age_prediction <- function(gse_id){
	mat = getting_geo_count_mat(gse_id)
	mat %<>% processing_geo_count_mat_for_age_prediction
	return(mat)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
getting_GEO_pdata <- function(gse_id, mat){
	lES = GEOquery::getGEO(gse_id, getGPL = FALSE)
	if(length(lES) > 1) stop('There is more than one data source. Manually 
		curation is needed') 
	ES = lES[[1]]
	raw_pdata = data.frame(Biobase::pData(lES[[1]]))
	return(raw_pdata)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
getting_GEO_ES <- function(mat, pdata){
	samples_mat = colnames(mat)
	samples_pda = rownames(pdata)
	if(!(length(samples_mat) == length(samples_pda))) stop('Sample length differs 
		between the count and annotation matrices')
	if(!all.equal(samples_mat, samples_pda)) stop('Sample ids differs between 
		the count and annotation matrices')
	ES = Biobase::ExpressionSet(mat, phenoData = Biobase::AnnotatedDataFrame(pdata))
	return(ES)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
getting_GEO_ES_for_age_model <- function(gse_id){
	mat = getting_geo_count_mat_for_age_prediction(gse_id)
	pdata = getting_GEO_pdata(gse_id)
	sel_samples = pdata$geo_accession %>% .[. %in% colnames(mat)]
	mat1 = mat[, sel_samples]
	pdata1 = pdata %>% .[match(sel_samples, .$geo_accession), ]
	len_pdata = nrow(pdata); len_mat = length(sel_samples)
	if(len_mat != len_pdata) print(paste(len_pdata - len_mat, 'samples', 
		'missing in the matrix'))
	ES = getting_GEO_ES(mat1, pdata1)
	return(ES)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
#' @examples
#' library(magrittr)
#' data(ES_GSE103938)
#' ES = ES_GSE103938
#' mat = Biobase::exprs(ES)
#' mat %<>% filtering_age_model_genes_and_rank_norm
#' ( v_age_scores = predicting_age_score(t(mat)) )
predicting_age_score <- function(mat, model_type = 'PASTA'){
	requireNamespace('glmnet', quietly = TRUE)
  data(cvfit_REG,    envir = environment())
  data(cvfit_PASTA,  envir = environment())
  data(beta_PASTA,   envir = environment())
  data(cvfit_C46,    envir = environment())
  data(beta_C46,     envir = environment())

	if(model_type == 'PASTA') cur_model = cvfit_PASTA
	if(model_type == 'REG')   cur_model = cvfit_REG
	if(model_type == 'CT46')  cur_model = cvfit_C46
	if(!model_type %in% c('PASTA', 'REG', 'CT46')) stop('Specify a valid model; 
		either PASTA, REG, or CT46')

	v_age_scores = stats::predict(cur_model, mat, s = 'lambda.min', 
		type = 'link')[,1] %>% as.numeric
	if(model_type == 'PASTA') v_age_scores = v_age_scores * beta_PASTA
	if(model_type == 'CT46')  v_age_scores = v_age_scores * beta_C46

	return(v_age_scores) 
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
adding_age_preds_to_pdata <- function(pdata, mat_t, REG = TRUE, PASTA = TRUE, 
	CT46 = FALSE){
	if(REG)   pdata[, REG   := predicting_age_score(mat_t, model_type = 'REG')]
	if(PASTA) pdata[, PASTA := predicting_age_score(mat_t, model_type = 'PASTA')]
	if(CT46)  pdata[, CT46  := predicting_age_score(mat_t, model_type = 'CT46')]
	colnames(pdata) %<>% gsub('\\.ch1', '', .)
	colnames(pdata) %<>% gsub('\\:ch1', '', .)
	colnames(pdata) %<>% gsub(' ', '_', .)
	colnames(pdata) %<>% gsub('\\.', '_', .)
	return(pdata)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
#' @examples
#' library(magrittr)
#' library(jsutil)
#' library(data.table)
#' ES = getting_GEO_ES_for_age_model('GSE103938') %T>% pdim # 8113 8
#' pdata = getting_pdata_with_age_scores(ES)
#' dcast(pdata, treated_with ~ vector, value.var = 'PASTA', fun.aggregate = mean)
#' dcast(pdata, treated_with ~ vector, value.var = 'REG', fun.aggregate = mean)
#' data(ES_GSE103938)
#' ES = ES_GSE103938
#' pdata = getting_pdata_with_age_scores(ES, filter_genes = TRUE, rank_norm = TRUE)
#' dcast(pdata, treated_with ~ vector, value.var = 'PASTA', fun.aggregate = mean)
getting_pdata_with_age_scores <- function(ES, filter_genes = FALSE, 
	rank_norm = FALSE, REG = TRUE, PASTA = TRUE, CT46 = FALSE){
	mat = Biobase::exprs(ES)
	if(filter_genes) mat %<>% filtering_age_model_genes_and_rank_norm
	if(rank_norm) mat %<>% applying_rank_normalization
	mat_t = t(mat)
	pdata = Biobase::pData(ES) %>% copy %>% setDT
	pdata %<>% adding_age_preds_to_pdata(mat_t, REG, PASTA, CT46)
	return(pdata)
}


#' Title of the function
#'
#' Description of the function.
#'
#' @export
filter_cell_types_in_seu_object <- function(seu, n_cell_min = 500, 
	dry_run = FALSE, verbose = FALSE){
	v_type_to_keep = seu$type %>% table %>% .[. >= n_cell_min] %>% names
	seu1 = seu[, seu$type %in% v_type_to_keep]
	if(verbose) print(table(seu1[[c('type', 'age')]]))
	if(!dry_run) return(seu1)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
making_pseudobulks_from_seurat <- function(seu, chunk_size = 1000, verbose = TRUE){
	
	requireNamespace('data.table', quietly = TRUE)
	requireNamespace('Seurat', quietly = TRUE)

	if(chunk_size == 1){

		seu_bulk = seu

	} else {

		pdata1 = seu[[c('type', 'age')]] %>% setDT
		pdata1[, id := 1:.N, c('type', 'age')]
		pdata1[, chunk := round(id / chunk_size) + 1]
		# pdata1[, .N / 100, type] ; pdata1[, max(chunk), type] # => the chunking works
		pdata1[, type_1 := paste0(type, '_', chunk)]
		seu$type_1 = pdata1$type_1

		n_types = length(unique(seu$type))
		if(n_types == 1){
			v_grouping_vars = c('type_1', 'age')
		} else {
			v_grouping_vars = c('type_1', 'age', 'type')
		}

		seu_bulk <- Seurat::AggregateExpression(seu, 
			group.by = v_grouping_vars, return.seurat = TRUE, verbose = FALSE)

		if(n_types == 1) seu_bulk$type = unique(seu$type)

	}

	seu_bulk$age %<>% as.numeric
	seu_bulk$chunk_size = chunk_size

  if(verbose) seu_bulk[[c('age', 'type')]] %>% table %>% t %>% print
	
	return(seu_bulk)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
predicting_age_from_pseudobulks <- function(seu_bulk, 
	REG = TRUE, PASTA = TRUE, CT46 = TRUE){

	requireNamespace('data.table', quietly = TRUE)
	requireNamespace('Seurat', quietly = TRUE)

  mat = Seurat::GetAssayData(seu_bulk, layer = 'data')
  if(class(mat) == 'dgCMatrix') mat %<>% as.matrix
	mat %<>% filtering_age_model_genes_and_rank_norm
	pdata = seu_bulk[[c('chunk_size', 'type', 'age')]] %>% setDT
	pdata %<>% adding_age_preds_to_pdata(t(mat), REG = REG, PASTA = PASTA, 
		CT46 = CT46)
  return(pdata)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
making_pseudobulks_and_predict_age <- function(seu, chunk_size = 1000, 
	verbose = TRUE, REG = TRUE, PASTA = TRUE, CT46 = TRUE){
	seu_bulk = making_pseudobulks_from_seurat(seu, chunk_size = chunk_size, 
		verbose = verbose)
	pdata = predicting_age_from_pseudobulks(seu_bulk, REG = REG, PASTA = PASTA, 
		CT46 = CT46)
	return(pdata)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
predicting_age_multiple_chunks <- function(seu, v_chunk_sizes = c(500, 1000), 
	REG = TRUE, PASTA = TRUE, CT46 = TRUE, verbose = TRUE){
	pdata_big = purrr::map(v_chunk_sizes, ~making_pseudobulks_and_predict_age(seu, .x, 
		REG = REG, PASTA = PASTA, CT46 = CT46)) %>% do.call(rbind, .)
	return(pdata_big)
}

#' Title of the function
#'
#' Description of the function.
#'
#' @export
#' @examples
#' library(magrittr)
#' library(jsutil)
#' data(seu_orozco_2020_retina_horizontal_cells)
#' seu = seu_orozco_2020_retina_horizontal_cells %T>% pdim # 57596 1875
#' rm(seu_orozco_2020_retina_horizontal_cells)
#' seu$age  = seu$development_stage %>% gsub('-year.*', '', .) %>% gsub('-', ' ', .)
#' seu$type = seu$cell_type %>% as.character
#' seu %>% filter_cell_types_in_seu_object(dry_run = TRUE, verbose = TRUE)
#' seu %<>% filter_cell_types_in_seu_object %T>% pdim # 57596 1875
#' seu_bulk = making_pseudobulks_from_seurat(seu)
#' ( pdata = predicting_age_from_pseudobulks(seu_bulk) )
#' ( pdata = making_pseudobulks_and_predict_age(seu) )
#' ( pdata_big = predicting_age_multiple_chunks(seu) )
#' get_cor_by_chunk_from_pdata_big(pdata_big)
get_cor_by_chunk_from_pdata_big <- function(pdata_big){
	dt_melt = data.table::melt(pdata_big, id.vars = c('chunk_size', 'type', 'age'), 
		variable.name = 'model_type', value.name = 'pred_age')
	dt_cor = dt_melt[, cor(age, pred_age), c('chunk_size', 'model_type')]
	dt_cor1 = data.table::dcast(dt_cor, chunk_size ~ model_type, value.var = 'V1')
	return(dt_cor1)
}

