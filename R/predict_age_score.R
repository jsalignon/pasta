
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

# predict_age_multiple_chunks <- function(v_chunk_sizes, seu, cvfit_REGR, 
# 	cvfit_AS40, cvfit_TC46){
# 	dt_age_cor1 = lapply(v_chunk_sizes, results_by_chunk_size, seu, cvfit_REGR, 
# 		cvfit_AS40, cvfit_TC46) %>% do.call(rbind, .)
# 	return(dt_age_cor1)
# }

