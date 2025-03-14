
# cd $PASTA_DIR

# R

# note: render each file from a fresh R session

source('doc/render_md.R')

rmarkdown::render('doc/Rmd/Liu_Polo_2020.Rmd', output_dir = 'doc/tutorials')
