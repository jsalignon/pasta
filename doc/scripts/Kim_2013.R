# Kim Aging Cell 2013 Analysis Tutorial using the Pasta package
# =============================================================
# This script demonstrates an analysis on dataset GSE41714 from Kim, Aging Cell, 2013.
# It covers:
#   - Data acquisition via GEOquery,
#   - Gene ID conversion using biomaRt,
#   - Processing expression data using Pasta functions,
#   - Adding age predictions,
#   - Reshaping data and computing correlations,
#   - Visualization.
#
# Ensure that the following packages are installed:
#   pasta, magrittr, data.table, ggplot2, gtools, GEOquery, biomaRt

# -------------------------------
# Load Required Libraries
# -------------------------------
library(pasta)
library(magrittr)
library(data.table)
library(ggplot2)
library(gtools)
library(GEOquery)
library(biomaRt)

# -------------------------------
# Define Output Files
# -------------------------------
file_Kim2013 <- 'output/ES_Kim2013.rds'
file_humanht12v4 <- 'output/dt_ids__humanht_12_v4.rds'

# Create output directory if needed
if (!dir.exists("output")) {
  dir.create("output", showWarnings = FALSE, recursive = TRUE)
}

# -------------------------------
# Data Acquisition
# -------------------------------
if (!file.exists(file_Kim2013)) {
  # Load the ExpressionSet for GSE41714 using GEOquery
  ES <- getGEO('GSE41714', GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  # Save the ExpressionSet for future use
  saveRDS(ES, file = file_Kim2013)
} else {
  ES <- readRDS(file_Kim2013)
}

# Extract the expression matrix and inspect phenotype data
mat <- exprs(ES)
print(dim(mat))
print(pData(ES))

# -------------------------------
# Gene ID Conversion using biomaRt
# -------------------------------
if (!file.exists(file_humanht12v4)) {
  
  # Load the expression matrix from the ExpressionSet
  mat <- exprs(ES)
  
  # Load Ensembl using the Pasta helper function
  ensembl <- get_mart_ensembl_human()
  check_biomart_attributes(ensembl, 'v4') 

  # If no results, try a different Ensembl version:
  ensembl_109 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 110)
  check_biomart_attributes(ensembl_109, 'v4')
  
  # Use an older version if needed:
  ensembl_108 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version = 105)
  check_biomart_attributes(ensembl_108, 'v4')
  # Expected: 'illumina_humanht_12_v4' should be listed.
  
  # Convert gene IDs using the provided function
  dt_ids <- converting_gene_ids_to_ensembl_gene_ids(mat, ensembl_108, 'illumina_humanht_12_v4')
  
  # Save conversion table for future use
  saveRDS(dt_ids, file = file_humanht12v4)
} else {
  dt_ids <- readRDS(file_humanht12v4)
}

# -------------------------------
# Processing Expression Data
# -------------------------------
# Reload ExpressionSet (if needed) and conversion table
ES <- readRDS(file_Kim2013)
dt_ids <- readRDS(file_humanht12v4)
mat <- exprs(ES)

# Rename rows in the expression matrix using the conversion table
# and filter genes using the Pasta pipeline
mat1 <- renaming_rows_in_mat_after_gene_id_conversion(mat, dt_ids) %T>% pdim
mat2 <- filtering_age_model_genes_and_rank_norm(mat1) %T>% pdim  # Expected dimensions: 8113 x 12

# Extract phenotype data; assume population doubling info is in column 31
pdata <- pData(ES) %>% setDT
pdata1 <- pdata[, 31] %>% set_colnames('population_doubling')

# -------------------------------
# Adding Metadata and Predicted Ages
# -------------------------------
# Define a table with doubling time and beta-gal measurements
dt_beta_gal <- data.table(
  doubling_time  = c(2, 2, 3, 5, 7, 10, 12, 14, 15, 20, 30, 40),
  nb_of_passages = c(16, 22, 54, 67, 75, 79, 81, 83, 85, 87, 88, 89),
  betagal_low    = c(0, 0, 0, 0, 12, 20, 42, 44, 45, 77, 79, 90),
  betagal_high   = c(0, 0, 0, 0, 0, 2, 7, 3, 8, 44, 45, 63),
  day            = c(25, 40, 110, 165, 210, 240, 265, 290, 320, 350, 405, 445)
)

# Check if the doubling times match
print(all.equal(as.numeric(pdata1$population_doubling[1:11]),
                dt_beta_gal$doubling_time[1:11]))  # should be TRUE

# Use dt_beta_gal as phenotype data and add age predictions to it
pdata1 <- copy(dt_beta_gal)
pdata1 %<>% adding_age_preds_to_pdata(t(mat2), CT46 = TRUE)
print(dim(pdata1))

# -------------------------------
# Reshaping to Long Format for Correlation Analysis
# -------------------------------
# Here we use the column names from dt_beta_gal as identifier variables.
v_ids <- names(dt_beta_gal)
pdata1_long <- melt(pdata1, id.vars = v_ids, variable.name = 'model_type')
pdata1_long <- pdata1_long[model_type == 'PASTA']
pdata1_long1 <- melt(pdata1_long, id.vars = c('model_type', 'value'),
                      variable.name = 'outcome', value.name = 'outcome_value')

# -------------------------------
# Computing Correlations
# -------------------------------
dt_cor <- pdata1_long1[, .(PCC = cor(value, outcome_value)), by = .(model_type, outcome)][order(-PCC)]
# Clean up outcome names
dt_cor[, outcome := outcome %>% gsub('_', ' ', .) %>% gsub('beta', 'B', .) %>% gsub('nb of ', '', .)]
dt_cor[, outcome := factor(outcome, levels = rev(as.character(dt_cor[model_type == 'PASTA']$outcome)))]
print(dt_cor)

# -------------------------------
# Visualization
# -------------------------------
p_cor_sen_day_2 <- ggplot(pdata1_long1[outcome == 'day'], aes(x = outcome_value, y = value)) + 
  xlab('Day') + 
  ylab('Age Score') + 
  geom_point(size = 1, alpha = 0.7) + 
  ggtitle('Senescence Timecourse') +
  theme_minimal()

print(p_cor_sen_day_2)
