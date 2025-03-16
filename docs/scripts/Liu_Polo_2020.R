# run_analysis.R
# This script runs the Pasta example analysis on dataset GSE149694
# from Liu & Polo (2020, Nature).
#
# Before running, ensure that the following packages are installed:
# - pasta
# - magrittr
# - data.table
# - ggplot2
# - gtools
#

# Load necessary libraries
library(pasta)
library(magrittr)
library(data.table)
library(ggplot2)
library(gtools)

# Create output directory if needed
if (!dir.exists("../output")) {
  dir.create("../output")
}

# =============================================================================
# Data Acquisition and Processing
# =============================================================================

# Get GEO ExpressionSet and Age Scores

# Download ExpressionSet for GSE149694; expect dimensions: 8113 x 8
file_LiuPolo2020 = '../output/ES_LiuPolo2020.rds'
if(!file.exists(file_LiuPolo2020)){
  # Download ExpressionSet for GSE149694
  ES <- getting_GEO_ES_for_age_model('GSE149694') %T>% pdim  # Dimensions: 8113 x 8
  saveRDS(ES, file_LiuPolo2020)
}

# Obtain phenotype data with age predictions and subset for plotting
pdata <- getting_pdata_with_age_scores(ES)
pdata1 <- pdata[, c('title', 'PASTA', 'REG', 'CT46')]


# =============================================================================
# Reshape and Clean the Data
# =============================================================================

# Reshape data and extract condition, time, and replicate information from the title.
dt <- melt(pdata1, id.vars = 'title', variable.name = 'model_type', value.name = 'age_score')
dt[, title := gsub('RNA-seq_', '', title)]
dt[, condition := gsub('-.*', '', title)]
dt[, time := gsub('.*-(.*)_rep.*', '\\1', title)]
dt[, time := factor(time, levels = mixedsort(unique(dt$time)))]
dt[, rep := gsub('.*rep', '', title)]
dt[, title := NULL]
dt[, time1 := as.integer(time)]


# =============================================================================
# Advanced Processing for Liu & Polo 2020 Data
# =============================================================================

# Load and Prepare Data
dt_pasta <- dt[model_type == 'PASTA']

# Separate fibroblast and non-fibroblast data
dt_fibro <- dt_pasta[condition == 'Fibroblast']
dt_liuPolo <- dt_pasta[condition != 'Fibroblast']

# For each non-fibroblast condition, append a copy of the fibroblast data with the condition overwritten
for(condition1 in unique(dt_liuPolo$condition)){
  dt_fibro1 <- copy(dt_fibro)[, condition := condition1]
  dt_liuPolo <- rbind(dt_liuPolo, dt_fibro1)
}

# Reformat time variables: remove letters and adjust specific values
dt_liuPolo[, time1 := time %>% 
  gsub('D', '', .) %>% 
  gsub('P3', '30', .) %>% 
  gsub('P10', '30', .) %>% 
  as.integer]

# Display unique condition and time mappings (this will print to the console)
print(unique(dt_liuPolo[, .(condition, time, time1)])[order(condition, time1)])


# =============================================================================
# Compute Correlations
# =============================================================================

# Calculate Pearson and Spearman correlations between age_score and time1, grouped by model_type and condition
dt_cor <- dt_liuPolo[, .(
  PCC = cor(age_score, time1, method = 'pearson'), 
  SCC = cor(age_score, time1, method = 'spearman')
), by = .(model_type, condition)][order(SCC)]

print(dt_cor)

# Display overall correlation values
cat("Pearson correlation between time1 and age_score:",
    cor(dt_liuPolo$time1, dt_liuPolo$age_score), "\n")

cat("Spearman correlation between time1 and age_score:",
    cor(dt_liuPolo$time1, dt_liuPolo$age_score, method = 'spearman'), "\n")


# =============================================================================
# Visualization
# =============================================================================

# Plot predicted age scores over time for the Liu & Polo dataset
p_LiuPolo2020 <- ggplot(dt_liuPolo, aes(x = time, y = age_score)) + 
  geom_point(size = 1, alpha = 0.7) +
  ggtitle('Reprogramming Timecourse') +
  theme_minimal()  # Replace with your custom theme if available

print(p_LiuPolo2020)


# =============================================================================
# Conclusion
# =============================================================================

# This script demonstrates how to use Pasta to process GEO datasets (GSE149694)
# and predict cellular age scores. The workflow includes data acquisition, reshaping,
# correlation analysis, and visualization of results.
