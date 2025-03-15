- [Introduction](#introduction)
- [Introduction](#introduction-1)
- [Data Acquisition, Age-Prediction, and
  Processing](#data-acquisition-age-prediction-and-processing)
  - [Loading Libraries and Dataset](#loading-libraries-and-dataset)
  - [Processing the Metadata](#processing-the-metadata)
  - [Filtering and Creating Pseudobulk
    Samples](#filtering-and-creating-pseudobulk-samples)
- [Results](#results)
  - [Correlation Analysis](#correlation-analysis)
  - [Visualization](#visualization)

# Introduction

# Introduction

This tutorial demonstrates an analysis using the **Pasta** package on
the dataset **L5 ET - MTG** from Gabitto, *Nature Neuroscience*, 2024.  
The study, *“Integrated multimodal cell atlas of Alzheimer’s disease”*,
is available at
[CellxGene](https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30)
and [Nature Neuroscience](https://doi.org/10.1038/s41593-024-01774-5).
In particular, the dataset “L5 ET - MTG: Seattle Alzheimer’s Disease
Atlas (SEA-AD)” of 2,590 cells was downloaded from CellxGene. These
cells corresponds to Layer 5 extratelencephalic cortical neurons in the
middle temporal gyrus.

# Data Acquisition, Age-Prediction, and Processing

In this section, we load the example Seurat object, process the
metadata, filter cell types, and create pseudobulk samples for age
prediction.

## Loading Libraries and Dataset

Loading libraries.

Loading data

``` r
data(seu_gabitto_2024_L5_ET_MTG_neuron)
seu <- seu_gabitto_2024_L5_ET_MTG_neuron %T>% pdim  # Dimensions: 36,412 genes x 2,590 cells
```

    ## [1] 36412  1392

``` r
rm(seu_gabitto_2024_L5_ET_MTG_neuron)
```

## Processing the Metadata

We adjust the metadata to extract age information and cell type
annotations.

``` r
# Process the metadata: remove extra characters from development_stage and convert cell type to character.
seu$age  <- seu$development_stage %>% gsub("-year.*", "", .) %>% 
  gsub("-", " ", .) %>% gsub('80 year old and over stage', '85', .)
seu$type <- seu$cell_type %>% as.character
```

Distribution of ages

``` r
seu@meta.data$age %>% table
```

    ## .
    ##  29  43  50  60  72  75  78  80  81  82  83  84  85  86  87  88  89 
    ##  54  16  77 161  11  26  36  74  22  72  91   7 583  46  48  41  27

Distribution of cell types

``` r
seu@meta.data$cell_type %>% table
```

    ## .
    ## L5 extratelencephalic projecting glutamatergic cortical neuron 
    ##                                                           1392

Preview cell type filtering (dry-run)

``` r
seu %>% filter_cell_types_in_seu_object(n_cell_min = 500, dry_run = TRUE, verbose = TRUE)
```

    ##                                                                 age
    ## type                                                              29  43  50
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  54  16  77
    ##                                                                 age
    ## type                                                              60  72  75
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 161  11  26
    ##                                                                 age
    ## type                                                              78  80  81
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  36  74  22
    ##                                                                 age
    ## type                                                              82  83  84
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  72  91   7
    ##                                                                 age
    ## type                                                              85  86  87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 583  46  48
    ##                                                                 age
    ## type                                                              88  89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  41  27

=\> This function allows to remove cell types without a given number of
cells.

## Filtering and Creating Pseudobulk Samples

We filter the Seurat object to retain only cell types with at least 500
cells, and then create pseudobulk samples for age prediction.

``` r
# Filter the Seurat object to keep only cell types with at least 500 cells
seu %<>% filter_cell_types_in_seu_object %T>% pdim
```

    ## [1] 36412  1392

``` r
# Predict age using a single chunk sizes
# dt_age_pred <- making_pseudobulks_and_predict_age(seu, chunk_size = 1000)
# => faster if needed

# Predict age using multiple pseudobulk chunk sizes
v_chunk_sizes <- 2^(0:10)
dt_age_pred <- predicting_age_multiple_chunks(seu, v_chunk_sizes)
```

    ##                                                                 age
    ## type                                                              29  43  50
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  54  16  77
    ##                                                                 age
    ## type                                                              60  72  75
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 161  11  26
    ##                                                                 age
    ## type                                                              78  80  81
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  36  74  22
    ##                                                                 age
    ## type                                                              82  83  84
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  72  91   7
    ##                                                                 age
    ## type                                                              85  86  87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 583  46  48
    ##                                                                 age
    ## type                                                              88  89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  41  27
    ##                                                                 age
    ## type                                                              29  43  50
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  28   9  39
    ##                                                                 age
    ## type                                                              60  72  75
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  81   7  14
    ##                                                                 age
    ## type                                                              78  80  81
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  19  38  12
    ##                                                                 age
    ## type                                                              82  83  84
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  37  47   5
    ##                                                                 age
    ## type                                                              85  86  87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 293  24  25
    ##                                                                 age
    ## type                                                              88  89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  21  15
    ##                                                                 age
    ## type                                                              29  43  50
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  15   5  20
    ##                                                                 age
    ## type                                                              60  72  75
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  41   4   7
    ##                                                                 age
    ## type                                                              78  80  81
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  10  19   7
    ##                                                                 age
    ## type                                                              82  83  84
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  19  24   3
    ##                                                                 age
    ## type                                                              85  86  87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 147  13  13
    ##                                                                 age
    ## type                                                              88  89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  11   8
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  8  3 11 21  2
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  4  5 10  4 10
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron 12  2 74  7  7
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  6  4
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  4  2  6 11  2
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  3  3  6  2  5
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  7  1 37  4  4
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  4  3
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  3  1  3  6  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  2  3  2  3
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  4  1 19  2  3
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  2
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  1  2  4  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  2  2  1  2
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  1 10  2  2
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  1
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  2  2  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  2  1  2
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  2  1  6  1  1
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  2  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  1  1
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  3  1  1
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  1  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  1  1
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  2  1  1
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1
    ##                                                                 age
    ## type                                                             29 43 50 60 72
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  1  1
    ##                                                                 age
    ## type                                                             75 78 80 81 82
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  1  1  1
    ##                                                                 age
    ## type                                                             83 84 85 86 87
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1  2  1  1
    ##                                                                 age
    ## type                                                             88 89
    ##   L5 extratelencephalic projecting glutamatergic cortical neuron  1  1

# Results

In this section, we compute correlations between the true age and
predicted age scores and visualize the results.

## Correlation Analysis

``` r
# Compute correlations by chunk size for different modeling strategies
dt_cor <- dt_age_pred[, .( 
  n_pseudobulks = .N,
  REG = cor(age, REG),
  PASTA = cor(age, PASTA),
  TC46 = cor(age, CT46)), by = chunk_size]
print(dt_cor)
```

    ##     chunk_size n_pseudobulks       REG     PASTA        TC46
    ##          <num>         <int>     <num>     <num>       <num>
    ##  1:          1          1392 0.0928153 0.1889565  0.01863358
    ##  2:          2           714 0.1205921 0.2478043  0.04883783
    ##  3:          4           366 0.2092726 0.3713957  0.03836583
    ##  4:          8           190 0.3391128 0.4838395  0.06944892
    ##  5:         16           104 0.4583045 0.6005741  0.04355135
    ##  6:         32            59 0.4755096 0.6534864 -0.01326929
    ##  7:         64            38 0.6432975 0.7552937  0.14246422
    ##  8:        128            27 0.6319542 0.7618285  0.15867524
    ##  9:        256            20 0.7581683 0.8369504 -0.02586498
    ## 10:        512            18 0.7734265 0.8399354  0.07893230
    ## 11:       1024            18 0.7765761 0.8408193  0.05759046

``` r
# Reshape the data into long format for plotting
cur_dt1 <- melt(dt_cor[, c(1, 3:5)], id.vars = "chunk_size", 
                variable.name = "Modeling_strategy", 
                value.name = "PCC")

# Optionally, add or adjust factor levels for the modeling strategies
model_levels <- c('REG', 'TC46', 'PASTA')
cur_dt1$Modeling_strategy <- factor(cur_dt1$Modeling_strategy, levels = model_levels)
print(cur_dt1)
```

    ##     chunk_size Modeling_strategy         PCC
    ##          <num>            <fctr>       <num>
    ##  1:          1               REG  0.09281530
    ##  2:          2               REG  0.12059205
    ##  3:          4               REG  0.20927264
    ##  4:          8               REG  0.33911276
    ##  5:         16               REG  0.45830445
    ##  6:         32               REG  0.47550962
    ##  7:         64               REG  0.64329745
    ##  8:        128               REG  0.63195424
    ##  9:        256               REG  0.75816830
    ## 10:        512               REG  0.77342651
    ## 11:       1024               REG  0.77657612
    ## 12:          1             PASTA  0.18895648
    ## 13:          2             PASTA  0.24780433
    ## 14:          4             PASTA  0.37139574
    ## 15:          8             PASTA  0.48383950
    ## 16:         16             PASTA  0.60057410
    ## 17:         32             PASTA  0.65348638
    ## 18:         64             PASTA  0.75529374
    ## 19:        128             PASTA  0.76182855
    ## 20:        256             PASTA  0.83695041
    ## 21:        512             PASTA  0.83993535
    ## 22:       1024             PASTA  0.84081929
    ## 23:          1              TC46  0.01863358
    ## 24:          2              TC46  0.04883783
    ## 25:          4              TC46  0.03836583
    ## 26:          8              TC46  0.06944892
    ## 27:         16              TC46  0.04355135
    ## 28:         32              TC46 -0.01326929
    ## 29:         64              TC46  0.14246422
    ## 30:        128              TC46  0.15867524
    ## 31:        256              TC46 -0.02586498
    ## 32:        512              TC46  0.07893230
    ## 33:       1024              TC46  0.05759046
    ##     chunk_size Modeling_strategy         PCC

## Visualization

We now plot the Pearson correlation coefficients (PCC) for the different
modeling strategies across chunk sizes.

``` r
# Plot PCC vs. log2(chunk_size) for different modeling strategies
p1 <- ggplot(cur_dt1, aes(x = log2(chunk_size), y = PCC, colour = Modeling_strategy)) +
  geom_point(size = 2) +
  geom_line() +
  scale_colour_manual(values = c('PASTA' = 'red2', 'REG' = 'dodgerblue', 'TC46' = 'forestgreen')) +
  ggtitle("Correlation between True Age and Predicted Age Scores") +
  xlab("log2(Chunk Size)") +
  ylab("Pearson Correlation Coefficient (PCC)") +
  theme_minimal()

print(p1)
```

![](gabitto_2024_files/figure-gfm/plot-correlations-1.png)<!-- -->
