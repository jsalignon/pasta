
# Tutorial for analyzing Seurat datasets with Pasta
# ========================================================
# This script demonstrates an analysis on the CellxGene collection 
# 1ca90a2d-2943-483d-b678-b809bf464c30, "L5 ET - MTG" from Gabitto, Nature 
# Neuroscience, 2024


# -------------------------------
# 1. Loading Libraries and Dataset
# -------------------------------
library(pasta)
library(jsutil)
library(magrittr)
library(data.table)
library(ggplot2)
library(Seurat)
# You can install any missing packages using `install.packages()` for CRAN 
# packages (`magrittr`, `data.table`, `ggplot2`, `Seurat`), or 
# `devtools::install_github()` for GitHub packages (`pasta`, `jsutil`).

# Create output directory if needed
if (!dir.exists("../output")) {
  dir.create("../output")
}


# -------------------------------
# 2. Downloading data
# -------------------------------

file_Gabitto2024 = '../output/seu_gabitto_2024.rds'
if(!file.exists(file_Gabitto2024)){
  download.file('https://datasets.cellxgene.cziscience.com/9d53f7bb-dc23-4c05-b2a6-4afa9a6e3be0.rds', 
    destfile = file_Gabitto2024)
}


# -------------------------------
# 3. Filtering data
# -------------------------------
# Filtering data to keep only cells from healthy donors with known age and 
# made with the 10X 3' method.

file_Gabitto2024_filtered = '../output/seu_gabitto_2024_filtered.rds'
if(!file.exists(file_Gabitto2024_filtered)){
  seu = readRDS(file_Gabitto2024)
  object.size(seu) # 572885136 bytes
  # Removing samples: removing 930 dementia patients
  seu = seu[, seu$disease == 'normal'] %T>% pncol
  object.size(seu) # 363434608 bytes
  # Removing 175 10x multiome samples,
  seu = seu[, seu$assay == '10x 3\' v3'] %T>% pncol
  object.size(seu) # 326000136 bytes
  # Removing 63 samples from donors or unknown age,
  seu = seu[, seu$development_stage != 'adult stage'] %T>% pncol
  object.size(seu) # 313914208 bytes
  saveRDS(seu, file_Gabitto2024_filtered)
}


# -------------------------------
# 4. Processing the Metadata
# -------------------------------
# Adjust metadata: extract age and convert cell type to character.
seu$age  <- seu$development_stage %>% gsub("-year.*", "", .) %>% 
  gsub("-", " ", .) %>% gsub('80 year old and over stage', '85', .)
seu$type <- seu$cell_type %>% as.character

# Distribution of ages
cat("Distribution of ages:\n")
print(table(seu@meta.data$age))

# Distribution of cell types
cat("Distribution of cell types:\n")
print(table(seu@meta.data$cell_type))

# Preview cell type filtering (dry-run)
cat("Preview cell type filtering (dry-run):\n")
seu %>% filter_cell_types_in_seu_object(n_cell_min = 500, dry_run = TRUE, verbose = TRUE)

# Keeping only cell types with at least 500 cells
seu %<>% filter_cell_types_in_seu_object %T>% pdim


# -------------------------------
# 5. Creating pseudobulks and prediting their age-effects
# -------------------------------

# Predict age using multiple pseudobulk chunk sizes.
v_chunk_sizes <- 2^(0:9)
dt_age_pred <- predicting_age_multiple_chunks(seu, v_chunk_sizes, verbose = F)


# -------------------------------
# 6. Correlation Analysis
# -------------------------------
# Compute correlations by chunk size for different modeling strategies.
dt_cor <- dt_age_pred[, .(
  n_pseudobulks = .N,
  REG = cor(age, REG),
  Pasta = cor(age, Pasta),
  TC46 = cor(age, CT46)
), by = chunk_size]
print(dt_cor)

# Reshape the correlation data into long format for plotting.
cur_dt1 <- melt(dt_cor[, c(1, 3:5)], id.vars = "chunk_size", 
                variable.name = "Modeling_strategy", 
                value.name = "PCC")

model_levels <- c("REG", "TC46", "Pasta")
cur_dt1$Modeling_strategy <- factor(cur_dt1$Modeling_strategy, levels = model_levels)


# -------------------------------
# 7. Visualization
# -------------------------------
# Plot Pearson correlation coefficients (PCC) for different modeling strategies.
p1 <- ggplot(cur_dt1, aes(x = log2(chunk_size), y = PCC, 
    colour = Modeling_strategy)) +
  geom_point(size = 2) +
  geom_line() +
  scale_colour_manual(values = c('Pasta' = 'red2', 
    'REG' = 'dodgerblue', 'TC46' = 'forestgreen')) +
  ggtitle("Correlation between True Age and Predicted Age Scores") +
  xlab("log2(Chunk Size)") +
  ylab("Pearson Correlation Coefficient (PCC)") +
  theme_minimal()

print(p1)
