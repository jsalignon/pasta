## Welcome

**Pasta** (Predicting **A**ge-**S**hift from **T**ranscriptomic **A**nalyses) is an R package designed to predict cellular age using transcriptomic data. The tool streamlines the process of preparing transcriptomic datasets, applying pre-trained models, and generating age predictions.

For detailed information on the underlying method, please refer to our [manuscript on bioRxiv](https://www.biorxiv.org/).

---

## Installation

Pasta is available on GitHub. You can install it directly using `devtools`:

```r
devtools::install_git("git@github.com:jsalignon/pasta.git", upgrade = "never")
```

Dependencies
It is recommended to also install the following packages:

magrittr – for piping and functional programming.
jsutil – for additional utility functions.
Install them as follows:

```r
install.packages("magrittr")
# jsutil is available on GitHub as well:
devtools::install_git("git@github.com:jsalignon/jsutil.git", upgrade = "never")
```

Other packages that should be installed to run all tutorials include data.table, ggplot2, gtools.


## Quick start

Predicting age from a GEO Series id.
 - Predicting age from a GEO Series id: [tutorial](doc/tutorials/Liu_Polo_2020.md), [script](doc/scripts/Liu_Polo_2020.R)  
In this case example we analyze a reprogramming bulk RNA-Seq timecourse dataset (LiuPolo2020).
 - Another...


## Example with a GEO dataset

Example of age prediction using only a GEO Series id. In this example, we analyze [GSE103938](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103938), with Title "Transcriptome analysis of OSKM- or RAS-induced senescent IMR90 fibroblasts treated with Rapamycin.".

```r
library(jsutil)
library(magrittr)
library(pasta)
library(data.table)

ES = getting_GEO_ES_for_age_model('GSE103938') %T>% pdim # 8113 8
pdata = get_pdata_with_age_scores(ES)
dcast(pdata, treated.with ~ vector, value.var = 'PASTA', fun.aggregate = mean)
dcast(pdata, treated.with ~ vector, value.var = 'REG', fun.aggregate = mean)
```

Another exapmple; GSE149694
# Liu..Polo 2020, Nature
# Reprogramming roadmap reveals route to human induced trophoblast stem cells 
# https://www.nature.com/articles/s41586-020-2734-6
# http://hrpi.ddnetbio.com/ ; GSE149694

```r
ES = getting_GEO_ES_for_age_model('GSE149694') %T>% pdim # 8113 8
pdata = getting_pdata_with_age_scores(ES)
pdata1 = pdata[, c('title', 'PASTA', 'REG', 'CT46')]

dt = melt(pdata1, id.vars = 'title', variable.name = 'model_type', 
	value.name = 'age_score')
dt[, title := gsub('RNA-seq_', '', title)]
dt[, condition := gsub('-.*', '', title)]
dt[, time := gsub('.*-(.*)_rep.*', '\\1', title)]
dt[, time := factor(time, levels = gtools::mixedsort(unique(dt$time)))]
dt[, rep := gsub('.*rep', '', title)]
dt[, title := NULL]
dt[, time1 := as.integer(time)]
```



## Example with a Seurat object
Here’s a brief example to get you started with Pasta using a Seurat object. In this example, we load sample data, filter cell types, aggregate into pseudobulk samples, and predict age scores.

```r
library(magrittr)
library(jsutil)
library(pasta)

# Load example dataset
data(seu_orozco_2020_retina_horizontal_cells)
seu <- seu_orozco_2020_retina_horizontal_cells %T>% pdim # dimensions: 57,596 genes x 1875 cells
rm(seu_orozco_2020_retina_horizontal_cells)

# Process the metadata
seu$age <- seu$development_stage %>% gsub("-year.*", "", .) %>% gsub("-", " ", .)
seu$type <- seu$cell_type %>% as.character

# Preview cell type filtering
seu %>% filter_cell_types_in_seu_object(dry_run = TRUE, verbose = TRUE)

# Filter the Seurat object to keep only cell types with at least 500 cells
seu %<>% filter_cell_types_in_seu_object %T>% pdim

# Create pseudobulks of 1000 cells each and predict age in one step
pdata <- making_pseudobulks_and_predict_age(seu, chunk_size = 1000)
print(pdata)

# Advanced: Predict age using multiple pseudobulk chunk sizes (500 and 1000)
pdata_big <- predicting_age_multiple_chunks(seu, v_chunk_sizes = c(500, 1000))
print(pdata_big)
```

## Documentation
For more detailed documentation and function reference, please check the reference manual in the repository.

If you have any questions or suggestions, feel free to open an issue or contact the maintainer.

## Citation
If you use Pasta in your research, please cite our manuscript:

[Manuscript Title], bioRxiv, [DOI link].

