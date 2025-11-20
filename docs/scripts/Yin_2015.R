
# Script for analyzing GEO mouse microarray datasets with Pasta
# =============================================================
# This script demonstrates an analysis on dataset GSE40156 from Yin, Methods Mol Biol, 2015.


# -------------------------------
# 1. Loading Libraries
# -------------------------------
library(magrittr)
library(data.table)
library(ggplot2)
library(gtools)
library(GEOquery)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(mouse4302.db)
library(pasta)
library(jsutil)
# You can install any missing packages using `install.packages()` for CRAN 
# packages (`magrittr`, `data.table`, `ggplot2`, `gtools`), 
# `BiocManager::install()` for bioconductor packages (`GEOquery`, `AnnotationDbi`, 
# `org.Mm.eg.db, mouse4302.db`), or `devtools::install_github()` for GitHub 
# packages (`pasta`, `jsutil`).

# Create output directory if needed
if (!dir.exists("docs/output")) {
  dir.create("docs/output")
}


# -------------------------------
# 2. Getting GEO ExpressionSet
# -------------------------------
file_Yin2015_mouse = 'docs/output/ES_Yin2015_mouse.rds'
if(!file.exists(file_Yin2015_mouse)){

  # Downloading a list of ExpressionSets from GEO
  lES <- getGEO('GSE40156', GSEMatrix = TRUE, getGPL = FALSE)

  # Keeping the Affymetrix Mouse Genome 430 2.0 Array (GPL1261) data
  ES <- lES[["GSE40156-GPL1261_series_matrix.txt.gz"]]

  saveRDS(ES, file_Yin2015_mouse)
} else {
  ES <- readRDS(file_Yin2015_mouse)
}
mat <- exprs(ES) %T>% pdim  # 45,101 probes, 42 samples


# -------------------------------
# 3. Converting gene ids using AnnotationDbi
# -------------------------------
file_dt_ids_affy_mouse430 <- "docs/output/dt_ids__affy_mouse430_2.rds"
if (!file.exists(file_dt_ids_affy_mouse430)) {

  dt_ids <- AnnotationDbi::select(
    mouse4302.db,
    keys    = rownames(mat),
    keytype = "PROBEID",
    columns = c("ENSEMBL")
  ) %>% setDT %T>% pnrow

  colnames(dt_ids) <- c("affy_mouse430_2", "ensembl_gene_id")

  # Keep only probes with an Ensembl identifier
  dt_ids <- dt_ids[!is.na(ensembl_gene_id), ] %T>% pnrow

  # Save conversion table for future use
  saveRDS(dt_ids, file = file_dt_ids_affy_mouse430)

} else {
  dt_ids <- readRDS(file_dt_ids_affy_mouse430)
}

# Rename rows of the matrix to Ensembl identifiers using helper from jsutil / pasta
mat <- renaming_rows_in_mat_after_gene_id_conversion(mat, dt_ids, 
    species = "mouse") %T>% pdim   # 18,353 genes, 42 samples


# -------------------------------
# 4. Preparing data for age-prediction
# -------------------------------
mat %<>% keeping_mouse_one2one_orthologs %T>% pdim # 3,010 genes, 42 samples
mat %<>% filtering_age_model_genes_and_rank_norm %T>% pdim # 8,113 genes, 42 samples


# -------------------------------
# 5. Creating metadata and predict ages
# -------------------------------
# Phenotype data from GEOquery ExpressionSet
pdata <- pData(ES) %>% setDT %>% copy
pdata %<>% .[, c("geo_accession", "age:ch1", "tissue:ch1", "strain:ch1")] 
colnames(pdata) <- c("geo_accession", "age", "tissue", "strain")

# Sanity check
all.equal(colnames(mat), pdata$geo_accession) # TRUE

# Adding Pasta age predictions
pdata %<>% adding_age_preds_to_pdata(t(mat)) %T>% pdim 

# Keeping only samples with non missing age
pdata <- pdata[!is.na(age)] %T>% pdim

# Age is counted in weeks; we convert it to integer
pdata[, age := gsub(" weeks", "", age) %>% as.integer]


# -------------------------------
# 6. Computing Correlations
# -------------------------------
dt_cor = pdata[, .(PCC = cor(age, Pasta)), c('tissue', 'strain')][order(-PCC)]
dt_cor[, strain := gsub('.*Apoe.*', 'Apoe', strain)]
dt_cor[, strain := gsub('.*C57B.*', 'WT', strain)]
dt_cor[, condition := paste0(tissue, ', ', strain)]
dt_cor[, condition := factor(condition, levels = dt_cor$condition %>% rev)]
print(dt_cor)


# -------------------------------
# 7. Visualization
# -------------------------------
p1 <- ggplot(dt_cor, aes(x = PCC, y = condition)) +
  geom_point(size = 5) +
  theme_minimal(base_size = 18) +
  ggtitle("Mouse microarray (GSE40156)") +
  ylab("Tissue, strain") +
  xlab("Pearson correlation coefficient")
print(p1)
