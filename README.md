# Pasta

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

Quick Start
Here’s a brief example to get you started with Pasta using a Seurat object. In this example, we load sample data, filter cell types, aggregate into pseudobulk samples, and predict age scores.

```r
library(magrittr)
library(jsutil)
library(pasta)

# Load example dataset
data(seu_orozco_2020_retina_horizontal_cells)
seu <- seu_orozco_2020_retina_horizontal_cells %T>% pdim # dimensions: 57596 x 1875
rm(seu_orozco_2020_retina_horizontal_cells)

# Process the metadata
seu$age <- seu$development_stage %>% gsub("-year.*", "", .) %>% gsub("-", " ", .)
seu$type <- seu$cell_type %>% as.character

# Preview cell type filtering
seu %>% filter_cell_types_in_seu_object(dry_run = TRUE, verbose = TRUE)

# Filter the Seurat object
seu %<>% filter_cell_types_in_seu_object %T>% pdim

# Create pseudobulk samples and predict age scores
seu_bulk <- making_pseudobulks_from_seurat(seu)
pdata <- predicting_age_from_pseudobulks(seu_bulk)
print(pdata)

# Alternative workflow: Create pseudobulks and predict age in one step
pdata <- making_pseudobulks_and_predict_age(seu)
print(pdata)

# Advanced: Predict age using multiple pseudobulk chunk sizes
pdata_big <- predicting_age_multiple_chunks(seu)
print(pdata_big)
```

Documentation
For more detailed documentation and function reference, please check the reference manual in the repository.

If you have any questions or suggestions, feel free to open an issue or contact the maintainer.

Citation
If you use Pasta in your research, please cite our manuscript:

[Manuscript Title], bioRxiv, [DOI link].

