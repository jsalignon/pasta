
utils::globalVariables(c('v_genes_model', 'v_new_names_geo_mat', 'cvfit_PASTA', 
	'cvfit_REG', 'coef_agediff'))


#' @export
#' @examples
#' check_if_url_is_valid(geo_url)
checking_if_url_is_valid <- function(url_in,t=2){
  con = url(url_in)
  check = suppressWarnings(try(open.connection(con, open = 'rt', timeout = t),
  	silent = T)[1])
  suppressWarnings(try(close.connection(con), silent = T))
  ifelse(is.null(check), TRUE, FALSE)
}

#' @export
getting_geo_url <- function(gse_id){
	urld    = 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts'
	geo_url = paste0(urld, '&acc=', gse_id, '&file=', gse_id, 
		'_raw_counts_GRCh38.p13_NCBI.tsv.gz')
	return(geo_url)
}

#' @export
#' @examples
#' getting_geo_valid_url('GSE121277')
getting_geo_valid_url <- function(gse_id){
	geo_url = getting_geo_url(gse_id)
	valid_url = checking_if_url_is_valid(geo_url)
  if(!valid_url) stop('Count matrix not available for this GSE id.')
	return(geo_url)
}

#' @export
#' @examples
#' mat = getting_geo_count_mat('GSE121276')
getting_geo_count_mat <- function(gse_id){
	urld <- 'https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts'
	path <- paste0(urld, '&acc=', gse_id, '&file=', gse_id, '_raw_counts_GRCh38.p13_NCBI.tsv.gz');
	geo_url = getting_geo_valid_url(gse_id)
	mat <- as.matrix(data.table::fread(geo_url, header = T, 
		colClasses = 'integer'), rownames = 1)
	return(mat)
}

#' @export
#' @examples
#' gse_id = 'GSE121276'
#' mat = getting_geo_count_mat(gse_id)
#' mat %<>% converting_geo_mat_gene_ids
converting_geo_mat_gene_ids <- function(mat){
	# Loading the PASTA model
  data(v_new_names_geo_mat, envir = environment())

	rownames(mat) = v_new_names_geo_mat
	mat = mat[grep('ENSG', v_new_names_geo_mat, value = T),]
	return(mat)

}

#' @export
#' @examples
#' # gse_id = 'GSE121276'
#' mat = getting_geo_count_mat(gse_id)
#' mat %<>% converting_geo_mat_gene_ids
#' mat %<>% filtering_geo_mat_genes_for_age_prediction(T) %>% dim
filtering_geo_mat_genes_for_age_prediction <- function(mat, 
	rank_normalization = T){
  data(v_genes_model, envir = environment())
	mat = mat[match(v_genes_model, rownames(mat)), ]
	median_value = median(c(mat), na.rm = T)
	mat[is.na(mat)] = median_value
	rownames(mat) = v_genes_model
	if(rank_normalization){
		mat %<>% applying_rank_normalization
	}
	return(mat)
}

#' @export
applying_rank_normalization <- function(mat, ties_method = 'average'){
	mat = apply(mat, 2, rank, ties.method = ties_method)
	return(mat)
}

#' @export
#' @examples
#' mat = getting_geo_count_mat_for_age_prediction('GSE121276')
getting_geo_count_mat_for_age_prediction <- function(gse_id, 
	rank_normalization = T){
	mat = getting_geo_count_mat(gse_id)
	mat %<>% converting_geo_mat_gene_ids
	mat %<>% filtering_geo_mat_genes_for_age_prediction(rank_normalization)
	return(mat)
}

#' @export
#' @examples
#' pdata = getting_GEO_pdata('GSE121276')
getting_GEO_pdata <- function(gse_id, mat){
	lES = GEOquery::getGEO(gse_id, getGPL = F)
	if(length(lES) > 1) stop('There is more than one data source. Manually 
		curation is needed') 
	ES = lES[[1]]
	raw_pdata = data.frame(Biobase::pData(lES[[1]]))
	return(raw_pdata)
}

#' @export
#' @examples
#' mat = getting_geo_count_mat_for_age_prediction('GSE121276')
#' pdata = getting_GEO_pdata('GSE121276')
#' ES = getting_GEO_ES(mat, pdata)
getting_GEO_ES <- function(mat, pdata){
	samples_mat = colnames(mat)
	samples_pda = rownames(pdata)
	if(!(length(samples_mat) == length(samples_pda))) stop('Sample length differs 
		between the count and annotation matrices')
	if(!all.equal(samples_mat, samples_pda)) stop('Sample ids differs between 
		the count and annotation matrices')
	ES = Biobase::ExpressionSet(mat, phenoData = AnnotatedDataFrame(pdata))
	return(ES)
}

#' @export
#' @importFrom magrittr %>% %<>% %T>%
#' @examples
#' gse_id = 'GSE103938'
#' ES = getting_GEO_ES_for_age_model('GSE121276') %T>% pdim # 8113 8
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

#' @export
#' @importFrom magrittr %>% %<>% %T>%
#' @examples
#' ES = getting_GEO_ES_for_age_model('GSE121276') %T>% pdim # 8113 8
#' v_age_scores = predicting_age_score(t(Biobase::exprs(ES))) )

predicting_age_score <- function(mat, model_type = 'PASTA'){
  data(cvfit_PASTA, envir = environment())
  data(cvfit_REG, envir = environment())
  data(coef_agediff, envir = environment())

	if(model_type == 'PASTA') cur_model = cvfit_PASTA
	if(model_type == 'REG')   cur_model = cvfit_REG
	if(!model_type %in% c('PASTA', 'REG')) stop('Specify a valid model; either 
		PASTA or REG')

	v_age_scores = predict(cur_model, mat, s = s_lambda, type = 'link')[,1] %>% 
		as.numeric
	if(model_type == 'PASTA') v_age_scores = v_age_scores * coef_agediff

	return(v_age_scores) 
}

#' @export
#' @importFrom data.table setDT := copy .SD dcast
#' @importFrom magrittr %>% %<>% %T>%
#' ES = getting_GEO_ES_for_age_model('GSE103938') %T>% pdim # 8113 8
#' pdata = get_pdata_with_age_scores(ES)
#' dcast(pdata, treated.with ~ vector, value.var = 'PASTA', fun.aggregate = mean)
#' dcast(pdata, treated.with ~ vector, value.var = 'REG', fun.aggregate = mean)
get_pdata_with_age_scores <- function(ES){
	mat_t = t(Biobase::exprs(ES))
	pdata = Biobase::pData(ES) %>% copy %>% setDT
	pdata[, REG   := predicting_age_score(mat_t, model_type = 'REG')]
	pdata[, PASTA := predicting_age_score(mat_t, model_type = 'PASTA')]
	colnames(pdata) %<>% gsub('\\.ch1', '', .)
	return(pdata)
}

