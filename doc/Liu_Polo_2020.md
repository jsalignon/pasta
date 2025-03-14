- [Introduction](#introduction)
- [Data Acquisition and Processing](#data-acquisition-and-processing)
  - [Get GEO ExpressionSet and Age
    Scores](#get-geo-expressionset-and-age-scores)
  - [Reshape and Clean the Data](#reshape-and-clean-the-data)
- [Advanced Processing for Liu & Polo 2020
  Data](#advanced-processing-for-liu-polo-2020-data)
  - [Load and Prepare Data](#load-and-prepare-data)
  - [Compute Correlations](#compute-correlations)
  - [Visualization](#visualization)
  - [Correlation Summary](#correlation-summary)
- [Conclusion](#conclusion)

# Introduction

This document provides an example analysis using Pasta on dataset
**GSE149694** from Liu & Polo (2020, *Nature*). The study, titled
*“Reprogramming roadmap reveals route to human induced trophoblast stem
cells”*, is available online at
[Nature](https://www.nature.com/articles/s41586-020-2734-6) and via
[HRPI](http://hrpi.ddnetbio.com/).

# Data Acquisition and Processing

In this section, we download the GEO ExpressionSet, obtain phenotype
data with age scores, and perform some data wrangling. The following
code executes the complete workflow using the Pasta functions.

## Get GEO ExpressionSet and Age Scores

``` r
# Download ExpressionSet for GSE149694
ES <- getting_GEO_ES_for_age_model('GSE149694') %T>% pdim  # Dimensions: 8113 x 8
```

    ## Features  Samples 
    ##     8113       32

``` r
# Obtain phenotype data with age predictions
pdata <- getting_pdata_with_age_scores(ES)
# Subset phenotype data for plotting
pdata1 <- pdata[, c('title', 'PASTA', 'REG', 'CT46')]
```

## Reshape and Clean the Data

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

# Advanced Processing for Liu & Polo 2020 Data

The following section demonstrates additional processing steps. We
adjust the fibroblast and non-fibroblast conditions and reassign time
values. Finally, we compute correlations and create plots.

## Load and Prepare Data

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
dt_liuPolo[, time0 := as.integer(time)]
dt_liuPolo[time == 'P10', time0 := 5]
# Display unique condition and time mappings
unique(dt_liuPolo[, c('condition', 'time', 'time0', 'time1')])[order(condition, time0)]
```

    ##     condition   time time0 time1
    ##        <char> <fctr> <int> <int>
    ##  1:     5iLAF     D3     1     3
    ##  2:     5iLAF     D7     2     7
    ##  3:     5iLAF    D13     3    13
    ##  4:     5iLAF    D21     4    21
    ##  5:      NHSM     D3     1     3
    ##  6:      NHSM     D7     2     7
    ##  7:      NHSM    D13     3    13
    ##  8:      NHSM    D21     4    21
    ##  9:    Primed     D3     1     3
    ## 10:    Primed     D7     2     7
    ## 11:    Primed    D13     3    13
    ## 12:    Primed    D21     4    21
    ## 13:    Primed     P3     5    30
    ## 14:    Primed    P10     5    30
    ## 15:      RSeT     D3     1     3
    ## 16:      RSeT     D7     2     7
    ## 17:      RSeT    D13     3    13
    ## 18:      RSeT    D21     4    21
    ## 19:   t2iLGoY     D3     1     3
    ## 20:   t2iLGoY     D7     2     7
    ## 21:   t2iLGoY    D13     3    13
    ## 22:   t2iLGoY    D21     4    21
    ##     condition   time time0 time1

## Compute Correlations

We calculate both Pearson and Spearman correlations between predicted
age scores and the true time values.

``` r
dt_cor <- dt_liuPolo[, .(
  PCC = cor(age_score, time0, method = 'pearson'), 
  SCC = cor(age_score, time0, method = 'spearman')
), by = c('model_type', 'condition')][order(SCC)]
print(dt_cor)
```

    ##    model_type condition        PCC        SCC
    ##        <fctr>    <char>      <num>      <num>
    ## 1:      PASTA      NHSM -0.9656991 -0.9759001
    ## 2:      PASTA     5iLAF -0.9540971 -0.9759001
    ## 3:      PASTA    Primed -0.9638784 -0.9320620
    ## 4:      PASTA   t2iLGoY -0.9209137 -0.8845194
    ## 5:      PASTA      RSeT -0.9125756 -0.8783101

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

![](doc/Liu_Polo_2020_files/figure-gfm/plot-timecourse-1.png)<!-- -->

## Correlation Summary

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

# Conclusion

This document demonstrates how to use Pasta to process GEO datasets
(here GSE149694) and predict cellular age scores. The workflow includes
data acquisition, reshaping, correlation analysis, and visualization of
results.

Feel free to modify the paths, themes, or parameters to suit your
analysis needs.
