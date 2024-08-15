
#' A function to return an object and printing its length.
#'
#' @param x A vector or list.
#' @return Returns the input and print its length.
#' @export
#' @examples
#' library(magrittr)
#' v1 = 1:10 %T>% plength
predict_age_shift <- function(mat1, pdata1, ctl_column, v_ctl_values, grp_column) {
    
    # Loading function, parameters and the model
    library(glmnet)
    library(jsutil)
    s_lambda          = 'lambda.min'
    predict_type_cmap = 'link'
    ties_method       = 'min'
    cvfit = readRDS('data/pasta.rds')
    final_genes_chr = rownames(cvfit$glmnet.fit$beta) %T>% plength # 8113

    # Matching genes with the input matrix
    ## we select the right genes and we rank normalize them
    pdata = pdata1
    # pdata = copy(pdata1)
    mat_trt = as.matrix(mat1)
    matched = match(final_genes_chr, rownames(mat_trt))
    print(paste('matched', length(which(!is.na(matched))), 'genes out of', 
        length(matched), 'genes'))
    mat_trt %<>% .[match(final_genes_chr, rownames(.)), ] %T>% pdim # 8113 216
    mat_trt %<>% rank_normalize_mat(ties_method) %>% t

    ## we compute controls for each condition
    which_ctl = which(pdata[[ctl_column]] %in% v_ctl_values)
    pdata_ctl = pdata[which_ctl]
    mat_ctl0  = mat_trt[which_ctl, ] %T>% pdim
    # pdata[pdata_ctl$id, ]$drug %>% unique # "DMSO"

    ## then we create a matrix the same size as treatment that contains the 
    #  appropriate control for each sample
    mat_ctl_plate = collapse_matrix_by_vector(mat = mat_ctl0, vec = pdata_ctl[[grp_column]])
    mat_ctl = mat_ctl_plate[pdata[[grp_column]], ] %T>% pdim # 

    ## finally, we can substract the right control to each sample
    mat_dif = mat_trt - mat_ctl

    ## finally we make age prediction and return the updated metadata table
    pdata[, age := predict(cvfit_final, mat_dif, s = s_lambda, 
        type = predict_type_cmap)[,1] %>% as.numeric]

    return(pdata)

}


# install_github('jsutil', username = 'jsalignon')
library(jsutil)

