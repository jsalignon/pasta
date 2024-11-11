

#' Test dataset GSE103938
#'
#' The raw counts ExpressionSet R object for GSE103938 obtained from GEO. Title of the dataset on GEO: "Transcriptome analysis of OSKM- or RAS-induced senescent IMR90 fibroblasts treated with Rapamycin.". 
#'
#' @format A data frame with 57232 genes and 21 samples
#' @source The object was obtained from the supplementary file GSE103938_GeneCounts.xlsx at the link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103938.
#'
"ES_GSE103938"


#' Te PASTA model
#'
#' The PASTA model is a glmnet model which was trained on transcriptomic samples of 21 datasets and multiple tissues to predict the age-shift between two pairs of samples. 
#'
#' @format A cv.glmnet created using the glmnet package.
#'
"cvfit_PASTA"


#' The regression model
#'
#' @format A cv.glmnet created using the glmnet package.
#'
"cvfit_REG"


#' Genes used to convert genes ids from GEO tables from entrez to ensembl.
#'
#' @format Vector of ensembl gene names.
#'
'v_new_names_geo_mat'


#' Genes used by the age-prediction models.
#'
#' @format Vector of ensembl gene names.
#'
'v_genes_model'


# ## creating the test dataset: GSE103938
# # -> 1 sample is missing: [1] "GSM2786649"... a control...
# # => needed to download the file "GSE103938_GeneCounts.xlsx" to finally get the 
# #    right count matrix...

# library(magrittr)
# library(Biobase)

# geo_id = 'GSE103938'

# # # loading the count table counts table from GEO
# # urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# # path <- paste0(urld, '&acc=', geo_id, '&file=', geo_id, '_raw_counts_GRCh38.p13_NCBI.tsv.gz');
# # mat <- as.matrix(data.table::fread(path, header = T, colClasses = "integer"), rownames = 1)

# mat1 = openxlsx::read.xlsx('GSE103938_GeneCounts.xlsx')
# mat = as.matrix(mat1[, -1]) %>% apply(2, as.integer) %>% set_rownames(mat1[,1])

# # # loading the gene annotation table
# # path  <- paste(urld, "file=Human.GRCh38.p13.annot.tsv.gz", sep = "&");
# # fdata <- data.table::fread(path, header = T) %>% as.data.frame
# # rownames(fdata) = fdata$GeneID

# # Loading the sample annotation table
# ES_tmp <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE)
# pdata = pData(ES_tmp[[1]])
# rownames(pdata) = pdata$title

# # all.equal(rownames(mat), rownames(fdata)) # T
# # mat = mat[rownames(fdata),]
# # all.equal(rownames(mat), rownames(fdata)) # T

# all.equal(colnames(mat), rownames(pdata)) # F

# ES_GSE103938 = ExpressionSet(mat, phenoData = AnnotatedDataFrame(pdata))

# usethis::use_data(ES_GSE103938)






# ## creating the test dataset: GSE266663
#  # -> 3 samples are missing!

# library(magrittr)

# # loading the count table counts table from GEO
# urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
# path <- paste(urld, "acc=GSE266663", "file=GSE266663_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep = "&");
# mat <- as.matrix(data.table::fread(path, header = T, colClasses = "integer"), rownames = 1)

# # loading the gene annotation table
# path  <- paste(urld, "file=Human.GRCh38.p13.annot.tsv.gz", sep = "&");
# fdata <- data.table::fread(path, header = T) %>% as.data.frame
# rownames(fdata) = fdata$GeneID

# # Loading the sample annotation table
# library(Biobase)
# ES_tmp <- GEOquery::getGEO('GSE266663', GSEMatrix = TRUE)
# pdata = pData(ES_tmp$GSE266663_series_matrix.txt.gz)

# all.equal(rownames(mat), rownames(fdata)) # T
# all.equal(colnames(mat), rownames(pdata)) # F
# mat = mat[, rownames(pdata)]
# all.equal(colnames(mat), rownames(pdata)) # F

# ES = ExpressionSet(mat, 
#     phenoData = AnnotatedDataFrame(pdata), 
#     fData = AnnotatedDataFrame(fdata))



# pData(ES_GSE266663)[,1]
# usethis::use_data(ES_GSE266663)


# tbl <- as.matrix(

# path <- paste(urld, "acc=GSE266663", "file=GSE255982_family.soft.gz", sep="&");
# data.table::fread(path, header=T)

#     GSE255982_family.soft.gz 4.0 Kb


#     Human.GRCh38.p13.annot.tsv.gz 2022-11-08 4.2 Mb

# # install_github('jsutil', username = 'jsalignon')
# library(jsutil)
