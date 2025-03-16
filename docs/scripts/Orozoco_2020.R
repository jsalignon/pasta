# Retina_Horizontal_Cells_Tutorial.R
# ====================================================
# Tutorial: Pasta Analysis on Human Retina Horizontal Cells
# Data: GSE135133 (Horizontal cells of human eye; 1,875 cells)
# Author: Jérôme Salignon
# Date: Sys.Date()
#
# This script demonstrates:
#   - Loading the example Seurat object,
#   - Processing metadata,
#   - Filtering cell types and creating pseudobulk samples,
#   - Predicting age scores using different pseudobulk chunk sizes,
#   - Computing correlations between true age and predicted age scores,
#   - Visualizing the results.
#
# Ensure that the following packages are installed:
#   pasta, magrittr, data.table, ggplot2, gtools, jsutil
# ====================================================

# -------------------------------
# 1. Load Required Libraries
# -------------------------------
library(magrittr)
library(jsutil)
library(pasta)
library(data.table)
library(ggplot2)
library(gtools)

# -------------------------------
# 2. Load Dataset
# -------------------------------
# Load the example dataset (horizontal cells of human retina)
data(seu_orozco_2020_retina_horizontal_cells)
seu <- seu_orozco_2020_retina_horizontal_cells %T>% pdim  # Dimensions: 57,596 genes x 1,875 cells
rm(seu_orozco_2020_retina_horizontal_cells)

data(seu_gabitto_2024_L5_ET_MTG_neuron)

seu = readRDS('data/9d53f7bb-dc23-4c05-b2a6-4afa9a6e3be0.rds')

# middle temporal gyrus
# L5 ET - MTG: Seattle Alzheimer's Disease Atlas (SEA-AD)
# L5 extratelencephalic projecting glutamatergic cortical neuron
seu_gabitto_2024_L5_ET_MTG_neuron = seu
save(seu_gabitto_2024_L5_ET_MTG_neuron, file = 'seu_gabitto_2024_L5_ET_MTG_neuron.rda')
seu_orozco_2020_retina_horizontal_cells.rda

seu = seu_gabitto_2024_L5_ET_MTG_neuron
seu = seu[, seu$disease == 'normal'] %T>% pdim
seu = seu[, seu$assay == '10x 3\' v3'] %T>% pdim
seu = seu[, seu$development_stage != 'adult stage'] %T>% pdim

# -------------------------------
# 3. Process Metadata
# -------------------------------
# Adjust metadata: extract 'age' and 'cell_type'
seu$age  <- seu$development_stage %>% gsub("-year.*", "", .) %>% 
  gsub("-", " ", .) %>% gsub('80 year old and over stage', '85', .)
seu$type <- seu$cell_type %>% as.character


# Display distribution of ages and cell types
cat("Distribution of ages:\n")
print(table(seu@meta.data$age))

cat("\nDistribution of cell types:\n")
print(table(seu@meta.data$cell_type))

# Preview cell type filtering (dry-run)
cat("\nPreview cell type filtering (dry-run):\n")
seu %>% filter_cell_types_in_seu_object(n_cell_min = 300, dry_run = TRUE, verbose = TRUE)


# -------------------------------
# 4. Filtering and Creating Pseudobulk Samples
# -------------------------------
# Filter the Seurat object: keep only cell types with at least 500 cells
seu %<>% filter_cell_types_in_seu_object(n_cell_min = 300) %T>% pdim # 57596 1834

# Create pseudobulk samples using a chunk size of 1000 cells
dt_age_pred <- making_pseudobulks_and_predict_age(seu, chunk_size = 1000)
cat("\nPseudobulk sample dimensions:\n")
print(dt_age_pred)

# Predict age using multiple pseudobulk chunk sizes
v_chunk_sizes <- 2^(0:10)
dt_age_pred <- predicting_age_multiple_chunks(seu, v_chunk_sizes)

dt_age_pred[, .(REG = cor(age, REG),
                              PASTA = cor(age, PASTA),
                              TC46 = cor(age, CT46)),
                          by = chunk_size]

dt_age_pred[, .N, chunk_size]

# -------------------------------
# 5. Correlation Analysis
# -------------------------------
# Compute correlations by chunk size for different modeling strategies
cur_dt_cor <- dt_age_pred[, .(REG = cor(age, REG),
                              PASTA = cor(age, PASTA),
                              TC46 = cor(age, CT46)),
                          by = chunk_size]

# Reshape the data into long format for plotting
cur_dt1 <- melt(cur_dt_cor, id.vars = "chunk_size",
                variable.name = "Modeling_strategy", value.name = "PCC")

# Optionally, set factor levels for the modeling strategies
model_levels <- c("REG", "TC46", "PASTA")
cur_dt1$Modeling_strategy <- factor(cur_dt1$Modeling_strategy, levels = model_levels)
cat("\nCorrelation data:\n")
print(cur_dt1)

# -------------------------------
# 6. Visualization
# -------------------------------
# Plot Pearson correlation coefficients (PCC) for the different modeling strategies
p1 <- ggplot(cur_dt1, aes(x = log2(chunk_size), y = PCC, colour = Modeling_strategy)) +
  geom_point(size = 2) +
  geom_line() +
  scale_colour_manual(values = c('PASTA' = 'red2',
                                 'REG' = 'dodgerblue',
                                 'TC46' = 'forestgreen')) +
  ggtitle("Correlation between True Age and Predicted Age Scores") +
  xlab("log2(Chunk Size)") +
  ylab("Pearson Correlation Coefficient (PCC)") +
  theme_minimal()

print(p1)
