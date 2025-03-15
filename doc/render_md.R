
render_Rmd_from_prefix <- function(prefix){
  rmarkdown::render(paste0('doc/Rmd/', prefix, '.Rmd'))

  file.rename(paste0('doc/Rmd/', prefix, '.md'), 
    paste0('doc/tutorials/', prefix, '.md'))

  output_folder = paste0('doc/tutorials/', prefix, '_files')
  if(dir.exists(output_folder)) unlink(output_folder, recursive = TRUE)
  file.rename(paste0('doc/Rmd/', prefix, '_files'), 
    paste0('doc/tutorials/', prefix, '_files'))
}

