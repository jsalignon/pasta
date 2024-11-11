
# #' A function to collapse replicates for the same control samples.
# #'
# #' @param mat_ctl0 Count matrix for the control samples.
# #' @param vec Character vector which indicates the label of each columns. Controls with the same label will be aggregated.
# #' @return Returns the input and print its length.
# #' @importFrom data.table setDT :=
# #' @export
# #' @examples
# #' library(magrittr)
# #' mat_ctl <- matrix(sample(c(rep(0, 10), 1:14)), nrow = 4, ncol = 6)
# #' colnames(mat_ctl) = paste0('ctl_', 1:6)
# #' rownames(mat_ctl) = paste0('gene_', 1:4)
# #' set.seed(1)
# #' mat_ctl0 = rank_normalize_mat(mat_ctl) %>% t
# #' vec = rep(c('ctl_A', 'ctl_B', 'ctl_C'), 2)
# #' collapse_matrix_by_vector(mat_ctl0, vec)
# collapse_matrix_by_vector <- function(mat_ctl0, vec){
#     dt = cbind(vec = vec, as.data.frame(mat_ctl0)) 
#     data.table::setDT(dt)
#     dt = dt[, lapply(.SD, mean), vec]
#     vec1 = as.character(dt$vec)
#     dt[, vec := NULL]
#     mat_ctl1 = as.matrix(dt)
#     # print(vec1)
#     rownames(mat_ctl1) = vec1
#     return(mat_ctl1)
# }


# #' A function to rank normalize the count matrix.
# #'
# #' @param mat1 Count matrix.
# #' @return Returns the input and print its length.
# #' @export
# #' @examples
# #' mat1 <- matrix(sample(c(rep(0, 10), 1:14)), nrow = 6, ncol = 4)
# #' set.seed(1)
# #' rank_normalize_mat(mat1)
# rank_normalize_mat <- function(mat1){
#     mat1 = apply(mat1, 2, rank, ties.method = 'min')
#     return(mat1)
# }


# utils::globalVariables(c('age_shift', 'cvfit_pasta'))

# #' A function to return an object and printing its length.
# #'
# #' @param mat1 Count matrix, with rownames being Ensembl gene ids.
# #' @param pdata1 Sample annotation table (i.e., phenotypic data).
# #' @param ctl_column The name of the column to use to select the control samples.
# #' @param v_ctl_values Samples with these values in the `ctl_column` will be used as controls.
# #' @param grp_column The name of the column to use to select groups of treatment and control samples to analyze together.
# #' @return Returns the input and print its length.
# #' @export
# #' @importFrom data.table setDT := copy .SD
# #' @importFrom magrittr %>% %<>% %T>%
# #' @importFrom stats predict
# #' @importFrom glmnet predict.glmnet
# #' @importFrom utils data
# #' @examples
# #' library(magrittr)
# #' library(jsutil)
# #' library(ggplot2)
# #' library(data.table)
# #' data(ES_GSE103938, envir = environment())
# #' mat1 = Biobase::exprs(ES_GSE103938)
# #' pdata1 = data.table(sample = colnames(mat1))
# #' pdata1[, type := gsub('_.*', '', sample)]
# #' pdata1[, rep  := gsub('.*_', '', sample)]
# #' ctl_column   = 'type'
# #' v_ctl_values = 'Prolif'
# #' grp_column   = 'rep'
# #' pdata = predict_age_shift(mat1, pdata1, ctl_column, v_ctl_values, grp_column)
# #' pdata[, type := factor(type, levels = 
# #'     c('Prolif', 'RAS', 'RAS1nM', 'RAS10nM', 'OSKM', 'OSKM1nM', 'OSKM10nM') %>% rev)]
# #' p1 = ggplot(pdata, aes(x = type, y = age_shift, fill = type)) + 
# #'     geom_bar(stat = 'identity') +
# #'     facet_grid(~rep) +
# #'     theme_classic() + 
# #'     coord_flip()
# #' # pdf('test.pdf'); print(p1); dev.off()
# #' # => Please observe that this dataset is not amazing as a control since the 
# #' #    OSKM label is misleading as the authors annotate these samples as OSKM-
# #' #    induced senescence. Usually, the age-shift effect of OSKM is much stronger.
# predict_age_shift <- function(mat1, pdata1, ctl_column, v_ctl_values, grp_column) {
    
#     # Loading the PASTA model
#     data(cvfit_pasta, envir = environment())
#     final_genes_chr = rownames(cvfit_pasta$glmnet.fit$beta) #%T>% plength # 8113

#     # Matching genes with the input matrix
#     ## we select the right genes and we rank normalize them
#     pdata = data.table::copy(pdata1)
#     mat_trt = as.matrix(mat1)
#     matched = match(final_genes_chr, rownames(mat_trt))
#     print(paste('matched', length(which(!is.na(matched))), 'genes out of', 
#         length(matched), 'genes'))
#     mat_trt = mat_trt[match(final_genes_chr, rownames(mat_trt)), ] #%T>% pdim # 8113 216
#     print(dim(mat_trt))
#     mat_trt %<>% rank_normalize_mat %>% t

#     ## Computing controls for each condition
#     which_ctl = which(pdata[[ctl_column]] %in% v_ctl_values)
#     pdata_ctl = pdata[which_ctl]
#     mat_ctl0  = mat_trt[which_ctl, ] #%T>% pdim
#     # pdata[pdata_ctl$id, ]$drug %>% unique # "DMSO"

#     ## Creating a matrix the same size as treatment that contains the
#     #  appropriate control for each sample
#     mat_ctl_plate = collapse_matrix_by_vector(mat_ctl0 = mat_ctl0, 
#         vec = pdata_ctl[[grp_column]])
#     mat_ctl = mat_ctl_plate[pdata[[grp_column]], ] # %T>% pdim # 

#     ## Substracting the right control to each sample
#     mat_dif = mat_trt - mat_ctl

#     ## Making age-shift predictions
#     pdata[, age_shift := predict(cvfit_pasta, mat_dif, s = 'lambda.min', 
#         type = 'link')[,1] %>% as.numeric]
#     pdata[, age_shift := age_shift * -3.21108]

#     # returning the updated metadata table
#     return(pdata)

# }

