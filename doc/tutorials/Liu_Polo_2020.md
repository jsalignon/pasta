- [Introduction](#introduction)
- [Data acquisition, age-prediction, and
  processing](#data-acquisition-age-prediction-and-processing)
  - [Loading all needed libraries](#loading-all-needed-libraries)
  - [Getting GEO ExpressionSet and predicted age
    scores](#getting-geo-expressionset-and-predicted-age-scores)
  - [Reshaping to long format and cleaning the
    metadata](#reshaping-to-long-format-and-cleaning-the-metadata)
  - [Adding initial time points](#adding-initial-time-points)
- [Results](#results)
  - [Computing correlations](#computing-correlations)
  - [Visualization](#visualization)

# Introduction

This document provides an example analysis using Pasta on dataset
**GSE149694** from Liu & Polo, *Nature*, 2020. The study, titled
*“Reprogramming roadmap reveals route to human induced trophoblast stem
cells”*, is available online at
[GEO]((https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149694)),
[Nature](https://www.nature.com/articles/s41586-020-2734-6), and on the
[HRPI website](http://hrpi.ddnetbio.com/).

# Data acquisition, age-prediction, and processing

In this section, we download the GEO ExpressionSet, obtain phenotype
data with age scores, and perform some data wrangling. The following
code executes the complete workflow using the Pasta functions.

## Loading all needed libraries

## Getting GEO ExpressionSet and predicted age scores

``` r
file_LiuPolo2020 = 'output/ES_LiuPolo2020.rds'
if(!file.exists(file_LiuPolo2020)){
  dir.create('output', showWarnings = F)
  # Download ExpressionSet for GSE149694
  ES <- getting_GEO_ES_for_age_model('GSE149694') %T>% pdim  # Dimensions: 8113 x 8
  saveRDS(ES, file_LiuPolo2020)
}

ES = readRDS(file_LiuPolo2020)
# Obtain phenotype data with age predictions
pdata <- getting_pdata_with_age_scores(ES)
# Subset phenotype data for plotting
pdata1 <- pdata[, c('title', 'PASTA', 'REG', 'CT46')]
```

## Reshaping to long format and cleaning the metadata

Here we reshape the data and extract information such as condition,
time, and replicate details from the title.

``` r
dt <- melt(pdata1, id.vars = 'title', variable.name = 'model_type', value.name = 'age_score')
dt[, title := gsub('RNA-seq_', '', title)]
dt[, condition := gsub('-.*', '', title)]
dt[, time := gsub('.*-(.*)_rep.*', '\\1', title)]
dt[, time := factor(time, levels = gtools::mixedsort(unique(dt$time)))]
dt[, rep := gsub('.*rep', '', title)]
dt[, title := NULL]
dt[, time1 := as.integer(time)]
```

## Adding initial time points

The initial time points (Fibroblast, day 3 and 7) is the same between
all conditions (‘NHSM’, ‘5iLAF’, ‘Primed’, ‘t2iLGoY’, ‘RSeT’). missing
for all samples. We add this time point to each condition so we can then
compute correlation.

``` r
dt_pasta <- dt[model_type == 'PASTA']

# Separate fibroblast and non-fibroblast data
dt_fibro <- dt_pasta[condition == 'Fibroblast']
dt_liuPolo <- dt_pasta[condition != 'Fibroblast']

# For each non-fibroblast condition, append a copy of the fibroblast data with the condition overwritten
for(condition1 in unique(dt_liuPolo$condition)){
  dt_fibro1 <- copy(dt_fibro)[, condition := condition1]
  dt_liuPolo <- rbind(dt_liuPolo, dt_fibro1)
}

# Reformat time variables
dt_liuPolo[, time1 := time %>% gsub('D', '', .) %>% gsub('P3', '30', .) %>% 
  gsub('P10', '30', .) %>% as.integer]
# Display unique condition and time mappings
unique(dt_liuPolo[, c('condition', 'time', 'time1')])[order(condition, time1)]
```

    ##     condition   time time1
    ##        <char> <fctr> <int>
    ##  1:     5iLAF     D3     3
    ##  2:     5iLAF     D7     7
    ##  3:     5iLAF    D13    13
    ##  4:     5iLAF    D21    21
    ##  5:      NHSM     D3     3
    ##  6:      NHSM     D7     7
    ##  7:      NHSM    D13    13
    ##  8:      NHSM    D21    21
    ##  9:    Primed     D3     3
    ## 10:    Primed     D7     7
    ## 11:    Primed    D13    13
    ## 12:    Primed    D21    21
    ## 13:    Primed     P3    30
    ## 14:    Primed    P10    30
    ## 15:      RSeT     D3     3
    ## 16:      RSeT     D7     7
    ## 17:      RSeT    D13    13
    ## 18:      RSeT    D21    21
    ## 19:   t2iLGoY     D3     3
    ## 20:   t2iLGoY     D7     7
    ## 21:   t2iLGoY    D13    13
    ## 22:   t2iLGoY    D21    21
    ##     condition   time time1

# Results

## Computing correlations

We calculate both Pearson and Spearman correlations between predicted
age scores and the true time values.

``` r
dt_cor <- dt_liuPolo[, .(
  PCC = cor(age_score, time1, method = 'pearson'), 
  SCC = cor(age_score, time1, method = 'spearman')
), by = c('model_type', 'condition')][order(SCC)]
print(dt_cor)
```

    ##    model_type condition        PCC        SCC
    ##        <fctr>    <char>      <num>      <num>
    ## 1:      PASTA      NHSM -0.9259211 -0.9759001
    ## 2:      PASTA     5iLAF -0.9017599 -0.9759001
    ## 3:      PASTA    Primed -0.9441355 -0.9320620
    ## 4:      PASTA   t2iLGoY -0.8548639 -0.8845194
    ## 5:      PASTA      RSeT -0.8461645 -0.8783101

We can also display overall correlation values.

``` r
cat("Pearson correlation between time1 and age_score:",
    cor(dt_liuPolo$time1, dt_liuPolo$age_score), "\n")
```

    ## Pearson correlation between time1 and age_score: -0.9045076

``` r
cat("Spearman correlation between time1 and age_score:",
    cor(dt_liuPolo$time1, dt_liuPolo$age_score, method = 'spearman'), "\n")
```

    ## Spearman correlation between time1 and age_score: -0.9541385

## Visualization

Finally, we visualize the correlation data and the timecourse of the age
predictions.

``` r
# Plot predicted age scores over time for the Liu & Polo dataset
p_LiuPolo2020 <- ggplot(dt_liuPolo, aes(x = time, y = age_score)) + 
  geom_point(size = 1, alpha = 0.7) +
  ggtitle('Reprogramming Timecourse') +
  theme_minimal()  # Replace with your custom theme if available

print(p_LiuPolo2020)
```

![](/raid/jersal/workspace/pasta/doc/tutorials/Liu_Polo_2020_files/figure-gfm/plot-timecourse-1.png)<!-- -->
