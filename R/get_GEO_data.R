
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
#' mat %<>% filtering_genes_and_rank_normalizing(T) %>% dim
filtering_genes_and_rank_normalizing <- function(mat, 
	rank_norm = T){
  data(v_genes_model, envir = environment())
	mat = mat[match(v_genes_model, rownames(mat)), ]
	median_value = median(c(mat), na.rm = T)
	mat[is.na(mat)] = median_value
	rownames(mat) = v_genes_model
	if(rank_norm){
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
	rank_norm = T){
	mat = getting_geo_count_mat(gse_id)
	mat %<>% converting_geo_mat_gene_ids
	mat %<>% filtering_genes_and_rank_normalizing(rank_norm)
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
	ES = Biobase::ExpressionSet(mat, phenoData = Biobase::AnnotatedDataFrame(pdata))
	return(ES)
}

#' @export
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
#' @examples
#' library(maggritr)
#' ES = getting_GEO_ES_for_age_model('GSE121276') %T>% pdim # 8113 8
#' v_age_scores = predicting_age_score(t(Biobase::exprs(ES))) )
predicting_age_score <- function(mat, model_type = 'PASTA'){
	library(glmnet)
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

	v_age_scores = predict(cur_model, mat, s = 'lambda.min', type = 'link')[,1] %>% 
		as.numeric
	if(model_type == 'PASTA') v_age_scores = v_age_scores * beta_PASTA
	if(model_type == 'CT46')  v_age_scores = v_age_scores * beta_C46

	return(v_age_scores) 
}

#' @export
add_age_preds_to_pdata <- function(pdata, mat_t, REG = T, PASTA = T, 
	CT46 = F){
	if(REG)   pdata[, REG   := predicting_age_score(mat_t, model_type = 'REG')]
	if(PASTA) pdata[, PASTA := predicting_age_score(mat_t, model_type = 'PASTA')]
	if(CT46)  pdata[, PASTA := predicting_age_score(mat_t, model_type = 'CT46')]
	colnames(pdata) %<>% gsub('\\.ch1', '', .)
	return(pdata)
}

#' @export
#' library(maggritr)
#' library(data.table)
#' ES = getting_GEO_ES_for_age_model('GSE103938') %T>% pdim # 8113 8
#' pdata = get_pdata_with_age_scores(ES)
#' dcast(pdata, treated.with ~ vector, value.var = 'PASTA', fun.aggregate = mean)
#' dcast(pdata, treated.with ~ vector, value.var = 'REG', fun.aggregate = mean)
get_pdata_with_age_scores <- function(ES, REG = T, PASTA = T, CT46 = F){
	mat_t = t(Biobase::exprs(ES))
	pdata = Biobase::pData(ES) %>% copy %>% setDT
	pdata %<>% add_age_preds_to_pdata(mat_t, REG, PASTA, CT46)
	return(pdata)
}

#' @export
#' data(seu_orozco_2020_retina_horizontal_cells)
#' seu = seu_orozco_2020_retina_horizontal_cells %T>% pdim # 57596 1875
#' rm(seu_orozco_2020_retina_horizontal_cells)
#' seu_bulk = make_pseudobulks_from_seurat(seu)
make_pseudobulks_from_seurat <- function(seu, chunk_size = 1000, verbose = T){
	pdata1 = seu[[c('type', 'true_age', 'age')]] %>% setDT %T>% pnrow # 2491
  pdata1[, id := 1:.N, c('type','true_age', 'age')]
  pdata1[, chunk := round(id / chunk_size) + 1]
  # pdata1[, .N / 100, type] ; pdata1[, max(chunk), type] # => the chunking works
  pdata1[, type_1 := paste0(type, '_', chunk)]

  seu$type_1 = pdata1$type_1
  seu_bulk <- Seurat::AggregateExpression(seu, 
  	group.by = c('type_1', 'true_age', 'age', 'type'), return.seurat = TRUE)

  if(verbose) seu_bulk[[c('age', 'type')]] %>% table %>% t

	return(seu_bulk)
}


#' @export
get_age_pred_for_pseudobulk <- function(seu, chunk_size1 = 1000, REG = T, 
	PASTA = T, CT46 = F){

	if(chunk_size1 == 1){
		seu_bulk = seu
	} else {
		seu_bulk = make_pseudobulks_from_seurat(seu, chunk_size1)
	}
  mat   = Seurat::GetAssayData(seu_bulk, layer = "data") %T>% pdim # 60579 479
  pdata = seu_bulk[[c('type', 'age')]] %>% setDT %T>% pnrow # 2491
  pdata$age %<>% as.numeric
  mat %<>% filtering_genes_and_rank_normalizing
  pdata %<>% add_age_preds_to_pdata(t(mat), REG, PASTA, CT46)
  pdata[, chunk_size := chunk_size1]
  return(pdata)
}


#' @export
predict_age_multiple_chunks <- function(v_chunk_sizes, seu){
	dt_age_cor1 = lapply(v_chunk_sizes, get_age_pred_for_pseudobulk, REG = T, 
		PASTA = T, CT46 = F) %>% do.call(rbind, .)
	return(dt_age_cor1)
}

#   dt_age_cor = data.table(
#     chunk_size = chunk_size,
#     n_uniq_val = nb_unique_values_per_sample,
#     n_zeros = nb_zeros,
#     REGR = cor(pdata1$age,  v_pred_REGR),
#     AS40 = cor(pdata1$age, -v_pred_AS40),
#     TC46 = cor(pdata1$age, -v_pred_TC46)
#     )


# seu = readRDS('../data/cellxgene/6934232e-4db4-423f-836d-bd730941aeba.rds') %T>% pdim #  60579 37121

# ##################
# # Lung aging atlas

# # Single-cell multiomic profiling of human lungs across age groups
# # LungMAP; 46,500 cells
# # https://cellxgene.cziscience.com/collections/625f6bf4-2f33-4942-962e-35243d284837
# # wget https://datasets.cellxgene.cziscience.com/6934232e-4db4-423f-836d-bd730941aeba.rds

# seu = readRDS('../data/cellxgene/6934232e-4db4-423f-836d-bd730941aeba.rds') %T>% pdim #  60579 37121

# v_true_age0 = seu$development_stage %>% gsub('-old.*', '', .) %>% gsub('-', ' ', .)
# table(v_true_age0) %>% t
#   #      3 year 31 year 31st week post fertilization human stage
#   # [1,]  19081   10393                                    17026
# v_true_age = gsub('31st.*', '31 weeks', v_true_age0)
# table(v_true_age) %>% t
# seu$true_age  = v_true_age
# v_age = v_true_age %>% gsub('31 weeks', '0.60', .) %>% gsub(' year', '', .) %>% as.numeric
# seu$age  = v_age
# # seu$true_age  = seu$development_stage %>% gsub('-old.*', '', .) %>% gsub('-', ' ', .)
# seu$type = seu$cell_type %>% as.character
# seu$true_age %>% table

# # filtering out conditions that have too few cells
# v_type_to_keep = seu$type %>% table %>% .[. > 500] %>% names
# seu = seu[, seu$type %in% v_type_to_keep] %T>% pdim # 20113   2491

# seu$type %>% table
# # alveolar macrophage              B cell       ciliated cell           club cell    endothelial cell     epithelial cell          fibroblast          macrophage            monocyte  myofibroblast cell
# #                7620                1034                 666                1291                2767                1219                7073                1148                2645                 558
# #         native cell natural killer cell            pericyte              T cell   type I pneumocyte  type II pneumocyte
# #                2206                1071                 803                3100                3618                7226

# seu[[c('true_age', 'type')]] %>% table %>% t
# #                      true_age
# # type                  3 year 31 weeks 31 year
# #   alveolar macrophage   4084      907    2629
# #   B cell                 514      317     203
# #   ciliated cell          341      207     118
# #   club cell              730      191     370
# #   endothelial cell      1594     1076      97
# #   epithelial cell        359      719     141
# #   fibroblast            2178     3906     989
# #   macrophage              81      855     212
# #   monocyte              1267      649     729
# #   myofibroblast cell     214      200     144
# #   native cell            762     1111     333
# #   natural killer cell    451      496     124
# #   pericyte               329      323     151
# #   T cell                1933      626     541
# #   type I pneumocyte     1064     1736     818
# #   type II pneumocyte    2802     1797    2627

# dt_age_lung_atlas = predict_age_multiple_chunks(
#   v_chunk_sizes, seu, cvfit_REGR, cvfit_AS40, cvfit_TC46)

# saveRDS(dt_age_lung_atlas, paste0(out_dir_PASTA, 'dt_age_lung_atlas.rds'))
# # dt_age_lung_atlas = readRDS(paste0(out_dir_PASTA, 'dt_age_lung_atlas.rds'))



# predict_age_by_chunk_size <- function(chunk_size, seu, cvfit_REGR, cvfit_AS40, cvfit_TC46){

#   pdata1 = seu[[c('type', 'true_age', 'age')]] %>% setDT %T>% pnrow # 2491
#   pdata1[, id := 1:.N, c('type','true_age', 'age')][, chunk := round(id / chunk_size) + 1]
#   # pdata1[, .N / 100, type] ; pdata1[, max(chunk), type] # => the chunking works
#   pdata1[, type_1 := paste0(type, '_', chunk)]

#   seu$type_1 = pdata1$type_1
#   seu_bulk <- Seurat::AggregateExpression(seu, 
#   	group.by = c('type_1', 'true_age', 'age', 'type'), return.seurat = TRUE)

#   seu_bulk[[c('age', 'type')]] %>% table %>% t

#   mat1   = Seurat::GetAssayData(seu_bulk, slot = "data") %T>% pdim # 60579 479
#   pdata1 = seu_bulk[[c('type', 'age')]] %>% setDT %T>% pnrow # 2491
#   pdata1$age %<>% as.numeric

#   v_genes_model = rownames(cvfit_AS40$glmnet.fit$beta)

#   matched = match(v_genes_model, rownames(mat1))

#   mat2 = mat1 %>% .[match(v_genes_model, rownames(.)), ] %T>% pdim # 8113 216
#   mat2[is.na(mat2)] = 0
#   mat3 = apply(mat2, 2, rank, ties.method = ties_method) %>% t
#   mat2_t = t(mat2)
#   mat2_t[is.na(mat2_t)] = 0

#   nb_unique_values_per_sample = mean(apply(mat2_t, 1, function(x) length(unique(x))))
#   nb_zeros = mean(apply(mat2_t, 1, function(x) length(which(x == 0))))

#   v_pred_REGR = predict(cvfit_REGR, mat3, s = s_lambda, type = 'link')[,1] %>% as.numeric
#   v_pred_AS40 = predict(cvfit_AS40, mat3, s = s_lambda, type = 'link')[,1] %>% as.numeric
#   v_pred_TC46 = predict(cvfit_TC46, mat3, s = s_lambda, type = 'link')[,1] %>% as.numeric

#   dt_age_cor = data.table(
#     chunk_size = chunk_size,
#     n_uniq_val = nb_unique_values_per_sample,
#     n_zeros = nb_zeros,
#     REGR = cor(pdata1$age,  v_pred_REGR),
#     AS40 = cor(pdata1$age, -v_pred_AS40),
#     TC46 = cor(pdata1$age, -v_pred_TC46)
#     )

#   return(dt_age_cor)

# }

