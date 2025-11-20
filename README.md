
## Welcome to Pasta!  

**Pasta** (P̲redicting a̲ge-s̲hift from t̲ranscriptomic a̲nalyses) is an R package designed to predict cellular age-effects from various kinds of human and mouse transcriptomic data, such as bulk and single-cell RNA-Seq, microarrays, and L1000 data. The tool streamlines the process of preparing transcriptomic datasets and making age-predictions using different models. For detailed information on the underlying method, please refer to our [preprint](https://www.biorxiv.org/content/10.1101/2025.06.04.657785v1).

<p align="center">
<img src="/docs/images/Figure_S9.png" width="800" />
</p>

Sections: [Installation](#installation), [Quick start](#quick-start), [Examples](#examples), [Citation](#citation), [Licence](#licence)


## Installation  

Pasta is available on GitHub. You can install it directly using `remotes`:
```r
install.packages("remotes")
remotes::install_github("jsalignon/pasta", dependencies = TRUE, upgrade = "never")
```


It is recommended to also install the following packages:
 - CRAN: magrittr, data.table, ggplot2, gtools, Seurat
 - Bioconductor: Biobase
 - GitHub: jsalignon/jsutil


## Quick start  

GEO human datasets with the label "Analyze with GEO2R" can be downloaded and preprocessed for age prediction in a single command. Here is an example for [GSE149694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149694) (an RNA-Seq reprogramming timecourse):
``` r
library(pasta)
library(data.table)
library(Biobase)
library(magrittr)
library(jsutil)
ES <- getting_GEO_ES_for_age_model('GSE149694') %T>% pdim  # Dimensions: 8113 x 8
pdata = getting_pdata_with_age_scores(ES, filter_genes = F, rank_norm = F)
```
The function `getting_GEO_ES_for_age_model` downloads and preprocesses the data (filter genes and do rank normalization).
The function `getting_pdata_with_age_scores` predict ages.

With an ExpressionSet object, one can preprocess the data and predict age-effects in a single line:
``` r
data(ES_GSE103938)
pdata = getting_pdata_with_age_scores(ES_GSE103938, filter_genes = T, rank_norm = T)
```

Here is how to predict age with a count matrix: 
``` r
mat = exprs(ES_GSE103938) %T>% pdim # 57,232 genes, 21 samples
pdata = pData(ES_GSE103938) %>% copy %>% setDT
mat %<>% filtering_age_model_genes_and_rank_norm %T>% pnrow # 8113 genes
pdata %<>% adding_age_preds_to_pdata(t(mat), REG = TRUE, Pasta = TRUE, CT46 = TRUE)
pdata[1:3, c('title', 'treated_with', 'vector', 'REG', 'Pasta', 'CT46')]
```

    ##        title  treated_with                 vector      REG     Pasta     CT46
    ##       <char>        <char>                 <char>    <num>     <num>    <num>
    ## 1:  Prolif_1          <NA>    empty vector (MSCV) 45.52050  7.179945 5.182232
    ## 2:    OSKM_1          <NA> OSKM expressing vector 49.22161 -4.438762 3.149072
    ## 3: OSKM1nM_1 1nM Rapamycin OSKM expressing vector 36.73217 -8.277933 4.138329

Age can be predicted using either REG, a regression model, CT46, a young vs old classifier with young/old cutoffs at <40 and >60 years, and Pasta, an age-shift model trained with pairs of samples from individuals of at least 40 years of age-difference.

``` r
print(dcast(pdata, treated_with ~ vector, value.var = 'Pasta', fun.aggregate = mean))
```

    ##      treated_with OSKM expressing vector RAS expressing vector empty vector (MSCV)
    ##            <char>                  <num>                 <num>               <num>
    ## 1:           <NA>              -3.816706              36.43337            4.762343
    ## 2: 10nM Rapamycin             -15.788035              33.73355                 NaN
    ## 3:  1nM Rapamycin              -6.103977              38.12827                 NaN
In this example, we can see the rejuvenating effect of OSKM and Rapamycin, and the aging effect of RAS-overexpression (i.e., Oncogene-induced senescence). Please note that Pasta predicts only relative age, so the predicted age-effects should mainly be compared between different samples of a given study.

Here are the results when using the regression model:
``` r
print(dcast(pdata, treated_with ~ vector, value.var = 'REG', fun.aggregate = mean))
```

    ##     treated_with OSKM expressing vector RAS expressing vector empty vector (MSCV)
    ##           <char>                  <num>                 <num>               <num>
    ##1:           <NA>               45.31800              56.65823            32.35902
    ##2: 10nM Rapamycin               23.29410              48.67683                 NaN
    ##3:  1nM Rapamycin               37.72791              58.28124                 NaN


## Examples  

Here are four examples with tutorials and scripts illustrating how to use Pasta: 

 - Example 1. Analyzing GEO human bulk RNA-Seq datasets with Pasta: [tutorial](docs/tutorials/Liu_Polo_2020.md), [script](docs/scripts/Liu_Polo_2020.R)  
In this example, we analyze a reprogramming bulk RNA-Seq timecourse dataset (Liu&Polo, 2020). 

 - Example 2. Analyzing GEO human microarray datasets with Pasta: [tutorial](docs/tutorials/Kim_2013.md), [script](docs/scripts/Kim_2013.R).
In this example, we analyze a senescence microarray timecourse dataset (Kim 2013). We use biomaRt to convert microarrays ids to ensembl gene ids.

 - Example 3. Analyzing human Seurat datasets with Pasta: [tutorial](docs/tutorials/gabitto_2024.md), [script](docs/scripts/gabitto_2024.R).
In this example, we analyze a single-cell transcriptomic dataset of cortical neurons (Gabitto 2024). We use different pseudobulk sizes to aggregate single-cell, predict age, and calculate Pearson correlation for each pseudobulk size.

 - Example 4. Analyzing GEO mouse microarray datasets with Pasta: [tutorial](docs/tutorials/Yin_2015.md), [script](docs/scripts/Yin_2015.R).


## Citation  

If you use Pasta in your research, please cite our manuscript:

Pasta, an age-shift transcriptomic clock, maps the chemical and genetic determinants of aging and rejuvenation, [Salignon et al, 2025, bioRxiv](https://www.biorxiv.org/content/10.1101/2025.06.04.657785v1).


## License
**Licence**: This source code is released under the MIT licence, included [here](LICENSE).
