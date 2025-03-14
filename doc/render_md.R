
# Define a post-processor function that modifies image paths
post_processor <- function(metadata, input, output, clean, verbose) {
  # Read the rendered output file (the .md file)
  text <- readLines(output, encoding = "UTF-8")
  # Replace image paths that begin with "Liu_Polo_2020_files/" by prefixing "doc/"
  # Adjust the pattern as needed; this uses a regular expression with perl = TRUE.
  text <- gsub("(!\\[[^\\]]*\\]\\()([^/].*?_files/)", "\\1doc/\\2", text, perl = TRUE)
  # Write the updated text back to the output file
  writeLines(text, output, useBytes = TRUE)
  return(output)
}

# Create a custom output format based on md_document with the post-processor attached
custom_md_document <- function(..., variant = "gfm") {
  fmt <- rmarkdown::md_document(..., variant = variant)
  fmt$post_processor <- post_processor
  fmt
}
