
# cd $PASTA_DIR

# R

# note: render each file from a fresh R session

source('doc/render_md.R')

rmarkdown::render('doc/Rmd/Liu_Polo_2020.Rmd', output_dir = 'doc/tutorials')
rmarkdown::render('doc/Rmd/Kim_2013.Rmd', output_dir = 'doc/tutorials')

rmarkdown::render('doc/Rmd/Liu_Polo_2020.Rmd')

  out <- rmarkdown::render(input)
out = rmarkdown::render('doc/Rmd/Liu_Polo_2020.Rmd')

render_Rmd_from_prefix <- function(prefix){
	rmarkdown::render(paste0('doc/Rmd/', prefix, '.Rmd'))

	file.rename(paste0('doc/Rmd/', prefix, '.md'), 
		paste0('doc/tutorials/', prefix, '.md'))

	output_folder = paste0('doc/tutorials/', prefix, '_files')
	if(dir.exists(output_folder)) unlink(output_folder, recursive = TRUE)
	file.rename(paste0('doc/Rmd/', prefix, '_files'), 
		paste0('doc/tutorials/', prefix, '_files'))
}


render_Rmd_from_prefix('Liu_Polo_2020')
render_Rmd_from_prefix('Kim_2013')

